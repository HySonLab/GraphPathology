"""Training-only target transform; no per-test-slide library-size statistics."""
import numpy as np

class CountTransform:
    """Standardize log1p(raw counts) using training rows only.

    This conservative first protocol does not divide by each test location's
    library size. It may retain depth effects, to be assessed at C1.
    """
    def fit(self, counts):
        z = self._log(counts)
        self.mean = z.mean(axis=0)
        self.scale = z.std(axis=0)
        self.scale[self.scale < 1e-8] = 1.0
        return self

    @staticmethod
    def _log(counts):
        c = np.asarray(counts)
        if not np.isfinite(c).all() or np.any(c < 0):
            raise ValueError('Counts must be finite and nonnegative')
        return np.log1p(c.astype(np.float64))

    def transform(self, counts):
        return ((self._log(counts) - self.mean) / self.scale).astype(np.float32)

    def observed_inputs(self, counts, mask):
        mask = np.asarray(mask)
        if mask.dtype != bool or mask.shape != (len(counts),):
            raise ValueError('Observation mask must be a boolean vector of length N')
        # Index BEFORE transforming: even NaN/Inf hidden labels cannot influence inputs.
        out = np.zeros(np.shape(counts), dtype=np.float32)
        if mask.any():
            out[mask] = self.transform(np.asarray(counts)[mask])
        return out, mask.copy()


def select_panels(training_counts, genes, sizes=(1, 50, 250, 1000)):
    """Training log-count variance ranking; intentionally not called formal HVGs."""
    variance = np.log1p(training_counts.astype(float)).var(axis=0)
    order = np.lexsort((np.asarray(genes), -variance))
    return {str(n): {'genes': np.asarray(genes)[order[:n]].tolist(),
                     'requested_size': n, 'actual_size': min(n, len(genes)),
                     'status': 'available' if n <= len(genes) else 'unavailable_panel_too_small'}
            for n in sizes}
