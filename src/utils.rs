use lazy_static::lazy_static;

lazy_static! {
    pub static ref PHRED_TO_ERROR: [f64; 128] = {
        let mut error_lookup: [f64; 128] = [1.0; 128];

        for i in 0..128 {
            if i >= 33 {
                error_lookup[i] = 10_f64.powf(-1.0 * ((i - 33) as f64) / 10.0);
            };
        }

        return error_lookup;
    };
}

#[inline]
pub fn mean_error_from_phred(phred: &[u8]) -> f64 {
    let error_sum: f64 = phred
        .iter()
        .map(|phred| {
            return PHRED_TO_ERROR[*phred as usize];
        })
        .sum::<f64>();

    let error_mean = error_sum / phred.len() as f64;

    error_mean
}
