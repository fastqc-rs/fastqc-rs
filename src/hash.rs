use ahash::AHasher;

/// A default Hasher for some faster hashing scheme (deterministic)
pub struct TxHasher(AHasher);

impl std::hash::Hasher for TxHasher {
    #[inline]
    fn finish(&self) -> u64 {
        self.0.finish()
    }
    #[inline]
    fn write(&mut self, bytes: &[u8]) {
        self.0.write(bytes)
    }
}

impl Default for TxHasher {
    fn default() -> Self {
        // just some random keys I came up with, need to be deterministic
        TxHasher(AHasher::new_with_keys(
            0x6c62_272e_07bb_0142,
            0x517c_c1b7_2722_0a95,
        ))
    }
}

/// A default BuildHasher using some faster hashing scheme (faster than SipHash)
pub type BuildHasher = std::hash::BuildHasherDefault<TxHasher>;

/// A default HashMap using some faster hashing scheme
pub type TxHashMap<K, V> = std::collections::HashMap<K, V, BuildHasher>;

/// A default HashSet using some faster hashing scheme
pub type TxHashSet<K> = std::collections::HashSet<K, BuildHasher>;
