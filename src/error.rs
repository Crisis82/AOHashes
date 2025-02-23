use dusk_safe::Error as SafeError;

/// Defines all possible error variants for SAFE
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Error {
    /// A call to during the lifetime of the [`safe::Sponge`] that doesn't fit
    /// the io-pattern.
    IOPatternViolation,

    /// An invalid io-pattern.
    InvalidIOPattern,

    /// The input doesn't yield enough input elements.
    TooFewInputElements,

    /// Failed to encrypt the message into the cipher with the provided secret
    /// and nonce.
    EncryptionFailed,

    /// Failed to decrypt the message from the cipher with the provided secret
    /// and nonce.
    DecryptionFailed,

    /// Invalid point on the jubjub-curve
    InvalidPoint,
}

impl From<SafeError> for Error {
    fn from(safe_error: SafeError) -> Self {
        match safe_error {
            SafeError::IOPatternViolation => Self::IOPatternViolation,
            SafeError::InvalidIOPattern => Self::InvalidIOPattern,
            SafeError::TooFewInputElements => Self::TooFewInputElements,
            SafeError::EncryptionFailed => Self::EncryptionFailed,
            SafeError::DecryptionFailed => Self::DecryptionFailed,
        }
    }
}
