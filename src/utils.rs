use lazy_static::lazy_static;

pub const PHRED_OFFSET: usize = 33;

lazy_static! {
    pub static ref PHRED_TO_ERROR: [f64; 128] = {
        let mut error_lookup: [f64; 128] = [1.0; 128];

        for (i, val) in error_lookup.iter_mut().enumerate() {
            if i >= PHRED_OFFSET {
                *val = 10_f64.powf(-((i - PHRED_OFFSET) as f64) / 10.0);
            };
        }

        error_lookup
    };
}

#[inline]
pub fn mean_error_from_phred(phred: &[u8]) -> f64 {
    let error_sum: f64 = phred
        .iter()
        .map(|phred| PHRED_TO_ERROR[*phred as usize])
        .sum::<f64>();

    error_sum / phred.len() as f64
}

#[cfg(test)]
mod tests {
    use super::*;
    use rstest::rstest;

    #[rstest]
    #[case(33, 1.0)]        // Phred 0 -> 100% error
    #[case(43, 0.1)]        // Phred 10 -> 10% error
    #[case(53, 0.01)]       // Phred 20 -> 1% error
    #[case(63, 0.001)]      // Phred 30 -> 0.1% error
    #[case(73, 0.0001)]     // Phred 40 -> 0.01% error
    fn test_phred_to_error(#[case] phred: usize, #[case] expected: f64) {
        let actual = PHRED_TO_ERROR[phred];
        assert!(
            (actual - expected).abs() < 1e-10,
            "PHRED_TO_ERROR[{}] = {}, expected {}",
            phred,
            actual,
            expected
        );
    }

    #[rstest]
    #[case(&[43, 43, 43, 43], 0.1)]           // All Phred 10 -> 10% mean error
    #[case(&[53, 53], 0.01)]                  // All Phred 20 -> 1% mean error
    #[case(&[43, 53], 0.055)]                 // Mixed: (0.1 + 0.01) / 2 = 0.055
    fn test_mean_error_from_phred(#[case] phred: &[u8], #[case] expected: f64) {
        let actual = mean_error_from_phred(phred);
        assert!(
            (actual - expected).abs() < 1e-10,
            "mean_error_from_phred({:?}) = {}, expected {}",
            phred,
            actual,
            expected
        );
    }

    #[test]
    fn test_phred_offset_values() {
        // Values below PHRED_OFFSET should be 1.0 (max error)
        for i in 0..PHRED_OFFSET {
            assert_eq!(
                PHRED_TO_ERROR[i], 1.0,
                "PHRED_TO_ERROR[{}] should be 1.0",
                i
            );
        }
    }
}
