use chrono::{DateTime, Local};
use itertools::Itertools;
use needletail::{parse_fastx_file, Sequence};
use rustc_hash::FxHashMap as HashMap;
use serde_json::json;
use serde_json::Value;
use std::error::Error;
use std::ffi::OsStr;
use std::fs::File;
use std::io;
use std::io::Write;
use std::path::Path;
use tera::{self, Context, Tera};

const BASES: [char; 5] = ['A', 'C', 'G', 'T', 'N'];
#[allow(unused)]
const A: usize = 0;
const C: usize = 1;
const G: usize = 2;
const T: usize = 3;
const N: usize = 4;

fn quartiles(hist: &[usize]) -> [f32; 5] {
    let sum = hist.iter().sum::<usize>();
    assert!(sum != 0);
    if sum == 1 {
        let value = hist.iter().enumerate().fold(
            0_usize,
            |acc, (value, &count)| if count > 0 { value } else { acc },
        ) as f32;
        return [value, value, value, value, value];
    }
    let mut ret = [0_f32; 5];
    // compute the quartiles
    for (i, quantile) in [0.25, 0.5, 0.75].iter().enumerate() {
        let rank = quantile * (sum - 1) as f64;
        let rank_ = rank.floor();
        let delta = rank - rank_;
        let n = rank_ as usize + 1;
        let mut acc = 0;
        let mut lo = None;
        for (hi, &count) in hist.iter().enumerate().filter(|(_, &count)| count > 0) {
            if acc == n && lo.is_some() {
                let lo = lo.unwrap() as f64;
                ret[i + 1] = (lo + (hi as f64 - lo) * delta) as f32;
                break;
            } else if acc + count > n {
                ret[i + 1] = hi as f32;
                break;
            }
            acc += count;
            lo = Some(hi);
        }
    }
    // compute lower, upper fences
    // TODO(lhepler): the UI reports these as min/max, which is incorrect.
    // If that's what we want, we can return those values.
    let iqr = ret[3] - ret[1];
    ret[0] = ret[1] - 1.5 * iqr;
    ret[4] = ret[3] + 1.5 * iqr;
    return ret;
}

pub(crate) fn process<P: AsRef<Path> + AsRef<OsStr>>(
    filename: P,
    k: u8,
    summary: Option<P>,
) -> Result<(), Box<dyn Error>> {
    let mut mean_read_qualities = HashMap::default();
    let mut base_count = HashMap::default();
    let mut base_quality_count = HashMap::default();
    let mut read_lengths = HashMap::default();
    let mut kmers = HashMap::default();
    let mut gc_content = vec![0_usize; 101];
    let mut read_count = 0_u64;
    let mut reader = parse_fastx_file(&filename).expect("Invalid path/file");
    let mut broken_read = false;

    // Gather data from every record
    while let Some(record) = reader.next() {
        if let Ok(seqrec) = record {
            read_count += 1;

            let sequence = seqrec.seq();
            let count_read_length = read_lengths.entry(sequence.len()).or_insert_with(|| 0_u64);
            *count_read_length += 1;

            let gc = sequence
                .iter()
                .filter(|&&b| b == b'G' || b == b'g' || b == b'C' || b == b'c')
                .count();
            let without_n = sequence.iter().filter(|&&b| b != b'N' || b != b'n').count();
            gc_content[gc * 100 / without_n] += 1;

            let mean_read_quality = if let Some(qualities) = seqrec.qual() {
                qualities.iter().map(|q| (q - 33) as u64).sum::<u64>() / qualities.len() as u64
            } else {
                0
            };

            let count_qual = mean_read_qualities
                .entry(mean_read_quality)
                .or_insert_with(|| 0);
            *count_qual += 1;

            if let Some(qualities) = seqrec.qual() {
                for (pos, &q) in qualities.iter().enumerate() {
                    let rec = base_quality_count
                        .entry(pos)
                        .or_insert_with(|| vec![0_usize; 94]);
                    rec[q as usize - 33] += 1;
                }
            }

            for (pos, base) in sequence.iter().enumerate() {
                let pos_count = base_count.entry(pos).or_insert_with(|| vec![0_u64; 5]);
                let i = match base {
                    b'A' | b'a' => 0,
                    b'C' | b'c' => 1,
                    b'G' | b'g' => 2,
                    b'T' | b't' => 3,
                    b'N' | b'n' => 4,
                    _ => panic!("Invalid base"),
                };
                pos_count[i] += 1;
            }

            let norm_seq = seqrec.normalize(false);
            let rc = norm_seq.reverse_complement();
            for (_, kmer, _) in norm_seq.canonical_kmers(k, &rc) {
                let count = kmers.entry(kmer.to_owned()).or_insert_with(|| 0_u64);
                *count += 1;
            }
        } else {
            broken_read = true;
        }
    }
    // Average gc content
    let avg_gc = {
        let (sum, len) = gc_content
            .iter()
            .enumerate()
            .fold((0_usize, 0_usize), |(s, l), (gc, n)| (s + gc * n, l + n));
        sum as f32 / len as f32
    };

    // Data for base per position
    let mut n_warn = "pass";
    let mut base_count_data = Vec::new();
    let mut base_count_percentage = HashMap::default();
    let mut gc_content_per_base = HashMap::default();
    let mut base_warning = "pass";
    for (position, bases) in base_count {
        let tmp_sum = bases.iter().sum::<u64>();
        let tmp_gc = (bases.get(G).unwrap() + bases.get(C).unwrap()) as f64
            / (tmp_sum - bases.get(N).unwrap()) as f64;
        gc_content_per_base.insert(position, tmp_gc);
        base_count_percentage.insert(
            position,
            bases
                .iter()
                .enumerate()
                .map(|(b, &count)| {
                    let (base, pct) = (BASES[b], count as f64 / tmp_sum as f64);
                    if base == 'N' && pct >= 20_f64 {
                        n_warn = "fail"
                    } else if base == 'N' && pct >= 5_f64 && n_warn != "fail" {
                        n_warn = "warn"
                    };
                    (base, pct)
                })
                .collect::<HashMap<char, f64>>(),
        );
        let gc_diff = i64::abs(*bases.get(G).unwrap() as i64 - *bases.get(C).unwrap() as i64)
            as f64
            / tmp_sum as f64;
        let tg_diff = i64::abs(*bases.get(T).unwrap() as i64 - *bases.get(G).unwrap() as i64)
            as f64
            / tmp_sum as f64;
        if gc_diff >= 20_f64 || tg_diff >= 20_f64 {
            base_warning = "fail"
        } else if (gc_diff >= 10_f64 || tg_diff >= 10_f64) && base_warning != "fail" {
            base_warning = "warn"
        }
        for (base, &count) in bases.iter().enumerate() {
            base_count_data.push(json!({"pos": position, "count": count, "base": BASES[base]}));
        }
    }

    let mut bpp_specs: Value =
        serde_json::from_str(include_str!("report/base_per_pos_specs.json"))?;
    bpp_specs["data"]["values"] = json!(base_count_data);

    // Data for read lengths
    let read_length_warn = if read_lengths.len() > 1 {
        "warn"
    } else {
        "pass"
    };
    let mut read_length_data = Vec::new();
    let mut read_length_sum = 0_u64;
    for (length, count) in &read_lengths {
        read_length_sum += *length as u64 * *count as u64;
        read_length_data.push(json!({"length": length, "count": count}));
    }
    let avg_read_length = read_length_sum / read_count as u64;

    let mut rle_specs: Value =
        serde_json::from_str(include_str!("report/read_lengths_specs.json"))?;
    rle_specs["data"]["values"] = json!(read_length_data);

    // Data for mean read qualities
    let mut mean_read_quality_data = Vec::new();
    for (quality, count) in mean_read_qualities {
        mean_read_quality_data.push(json!({"score": quality, "count": count}))
    }

    let mut sqc_specs: Value =
        serde_json::from_str(include_str!("report/sequence_quality_score_specs.json"))?;
    sqc_specs["data"]["values"] = json!(mean_read_quality_data);

    // Data for kmer quantities
    let mut kmer_data = Vec::new();
    for (kmer, count) in &kmers {
        kmer_data.push(json!({"k_mer": std::str::from_utf8(kmer).unwrap(), "count": count}))
    }

    let mut overly_represented = Vec::new();
    let mut overly_represented_warn = "pass";
    let total_kmers = kmers.iter().map(|(_, count)| count).sum::<u64>();
    for (km, occ) in kmers
        .iter()
        .sorted_by(|(_, a), (_, b)| Ord::cmp(&b, &a))
        .take(5)
    {
        let percentage = *occ as f64 / total_kmers as f64;
        if percentage >= 1_f64 {
            overly_represented_warn = "fail";
            overly_represented.push(json!({"k_mer": std::str::from_utf8(km).unwrap(), "count": occ, "pct": percentage, "or": "Yes"}));
        } else if percentage >= 0.2_f64 {
            if overly_represented_warn != "fail" {
                overly_represented_warn = "warn";
            }
            overly_represented.push(json!({"k_mer": std::str::from_utf8(km).unwrap(), "count": occ, "pct": percentage, "or": "No"}));
        };
    }

    let mut counter_specs: Value = serde_json::from_str(include_str!("report/counter_specs.json"))?;
    counter_specs["data"]["values"] = json!(kmer_data);

    // Data for GC content
    let mut gc_data = Vec::new();
    // TODO(lhepler): I had to add this filter to get the same histogram in the output report
    // but I think it may be more correct to leave it out. Arguable, though
    for (perc, &count) in gc_content
        .iter()
        .enumerate()
        .filter(|(_, &count)| count > 0)
    {
        gc_data.push(json!({"gc_pct": perc, "count": count, "type": "gc"}))
    }

    let mut gc_specs: Value = serde_json::from_str(include_str!("report/gc_content_specs.json"))?;
    gc_specs["data"]["values"] = json!(gc_data);

    // Data for base quality per position
    let mut base_quality_warn = "pass";
    let mut base_per_pos_data = Vec::new();
    for (position, qualities) in base_quality_count {
        let (sum, len) = qualities
            .iter()
            .enumerate()
            .fold((0_usize, 0_usize), |(s, l), (q, c)| (s + q * c, l + c));
        let avg = sum as f64 / len as f64;
        let values = quartiles(&qualities);
        if values.get(2).unwrap() <= &20_f32 {
            base_quality_warn = "fail"
        } else if values.get(2).unwrap() <= &25_f32 && base_quality_warn != "fail" {
            base_quality_warn = "warn"
        }
        base_per_pos_data.push(json!({
        "pos": position,
        "average": avg,
        "upper": values.get(4).unwrap(),
        "lower": values.get(0).unwrap(),
        "q1": values.get(1).unwrap(),
        "q3": values.get(3).unwrap(),
        "median":values.get(2).unwrap(),
        }));
    }

    let mut qpp_specs: Value =
        serde_json::from_str(include_str!("report/quality_per_pos_specs.json"))?;
    qpp_specs["data"]["values"] = json!(base_per_pos_data);

    let plots = json!({
        "k-mer quantities": {"short": "count", "specs": counter_specs.to_string()},
        "gc content": {"short": "gc", "specs": gc_specs.to_string()},
        "base sequence quality": {"short": "base", "specs": qpp_specs.to_string()},
        "sequence quality score": {"short": "qual", "specs": sqc_specs.to_string()},
        "base sequence content": {"short": "cont", "specs": bpp_specs.to_string()},
        "read lengths": {"short": "rlen", "specs": rle_specs.to_string()},
    });

    let file = Path::new(&filename).file_name().unwrap().to_str().unwrap();
    let meta = json!({
        "file name": {"name": "file name", "value": file},
        "canonical": {"name": "canonical", "value": "True"},
        "k": {"name": "k", "value": k},
        "total reads": {"name": "total reads", "value": read_count},
        "average GC content": {"name": "average GC content", "value": avg_gc},
        "average read length": {"name": "average read length", "value": avg_read_length},
    });

    let mut templates = Tera::default();
    templates.register_filter("embed_source", embed_source);
    templates.add_raw_template("report.html.tera", include_str!("report/report.html.tera"))?;
    let mut context = Context::new();
    context.insert("plots", &plots);
    context.insert("meta", &meta);
    let local: DateTime<Local> = Local::now();
    context.insert("time", &local.format("%a %b %e %T %Y").to_string());
    context.insert("version", &env!("CARGO_PKG_VERSION"));
    context.insert("invalid_reads", &broken_read);
    let html = templates.render("report.html.tera", &context)?;
    io::stdout().write_all(html.as_bytes())?;

    if let Some(path) = summary {
        let output_path = Path::new(&path);
        templates.add_raw_template(
            "fastqc_summary.txt.tera",
            include_str!("report/fastqc_summary.txt.tera"),
        )?;
        context.insert("filename", &file);
        context.insert("reads", &read_count);
        context.insert("avg_read_length", &avg_read_length);
        context.insert("avg_gc", &avg_gc);
        context.insert("bpp_data", &base_per_pos_data);
        context.insert("mean_read_quality_data", &mean_read_quality_data);
        context.insert("base_count", &base_count_percentage);
        context.insert("base_warning", &base_warning);
        context.insert("read_lengths", &read_lengths);
        context.insert("read_length_warn", &read_length_warn);
        context.insert("gc_data", &gc_data);
        context.insert("gc_per_base", &gc_content_per_base);
        context.insert("overly_represented", &overly_represented);
        context.insert("overly_represented_warn", &overly_represented_warn);
        context.insert("n_warn", &n_warn);
        context.insert("base_quality_warn", &base_quality_warn);
        let txt = templates.render("fastqc_summary.txt.tera", &context)?;
        let mut file = File::create(output_path.join("fastqc_data.txt"))?;
        file.write_all(txt.as_bytes())?;
    }
    Ok(())
}

fn embed_source(
    value: &tera::Value,
    _: &std::collections::HashMap<String, tera::Value>,
) -> tera::Result<tera::Value> {
    let url = tera::try_get_value!("upper", "value", String, value);
    let source = reqwest::get(&url).unwrap().text().unwrap();
    Ok(tera::to_value(source).unwrap())
}

#[cfg(test)]
mod test {
    use super::quartiles;
    #[test]
    fn test_quartiles1() {
        let v1 = [-49.5, 24.75, 49.5, 74.25, 148.5];
        let v2 = quartiles(&vec![1; 100]);
        assert!(v1 == v2);
    }
    #[test]
    fn test_quartiles2() {
        let v1 = [6.25, 25.0, 25.0, 37.5, 56.25];
        let mut h = [0_usize; 76];
        h[25] = 75;
        h[75] = 25;
        let v2 = quartiles(&h);
        assert!(v1 == v2);
    }
    #[test]
    fn test_quartiles3() {
        let v1 = [43.75, 62.5, 75.0, 75.0, 93.75];
        let mut h = [0_usize; 76];
        h[25] = 25;
        h[75] = 75;
        let v2 = quartiles(&h);
        assert!(v1 == v2);
    }
}
