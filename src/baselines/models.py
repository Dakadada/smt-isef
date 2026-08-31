"""Deterministic weighted ridge regression with an effective-DOF budget."""
import math


def _inverse(a):
    n = len(a)
    aug = [row[:] + [float(i == j) for j in range(n)] for i, row in enumerate(a)]
    for col in range(n):
        pivot = max(range(col, n), key=lambda r: abs(aug[r][col]))
        aug[col], aug[pivot] = aug[pivot], aug[col]
        if abs(aug[col][col]) < 1e-12:
            raise ValueError("singular ridge design")
        scale = aug[col][col]
        aug[col] = [v / scale for v in aug[col]]
        for row in range(n):
            if row != col:
                scale = aug[row][col]
                aug[row] = [v - scale * x for v, x in zip(aug[row], aug[col])]
    return [row[n:] for row in aug]


def _matvec(a, b):
    return [sum(x * y for x, y in zip(row, b)) for row in a]


class CapacityMatchedLinear:
    """Weighted full-representation ridge capped by true hat-matrix EDF."""
    def __init__(self, max_effective_df=8.0, ridge=1.0):
        self.max_effective_df, self.ridge = float(max_effective_df), float(ridge)
        if self.max_effective_df < 1 or self.ridge <= 0:
            raise ValueError("EDF budget must be >=1 and ridge must be positive")

    def fit(self, rows, y, weights=None):
        if not rows or len(rows) != len(y):
            raise ValueError("nonempty aligned rows and labels required")
        weights = list(weights or [1.0] * len(rows))
        if len(weights) != len(rows) or any(w <= 0 for w in weights):
            raise ValueError("positive aligned weights required")
        self.names = sorted({k for row in rows for k in row})
        total_w = sum(weights)
        self.means = {k: sum(w * row.get(k, 0.0) for row, w in zip(rows, weights)) / total_w for k in self.names}
        self.scales = {}
        for key in self.names:
            variance = sum(w * (row.get(key, 0.0) - self.means[key]) ** 2 for row, w in zip(rows, weights)) / total_w
            self.scales[key] = math.sqrt(variance) or 1.0
        x = [[1.0] + [(row.get(k, 0.0) - self.means[k]) / self.scales[k] for k in self.names] for row in rows]
        p = len(x[0])
        gram = [[sum(w * row[i] * row[j] for row, w in zip(x, weights)) for j in range(p)] for i in range(p)]
        rhs = [sum(w * row[i] * target for row, target, w in zip(x, y, weights)) for i in range(p)]
        penalized = [row[:] for row in gram]
        for i in range(1, p):
            penalized[i][i] += self.ridge
        inverse = _inverse(penalized)
        self.coef = _matvec(inverse, rhs)
        self.effective_degrees_of_freedom = sum(inverse[i][j] * gram[j][i] for i in range(p) for j in range(p))
        if self.effective_degrees_of_freedom > self.max_effective_df + 1e-8:
            raise ValueError("selected ridge exceeds the pre-registered effective-DOF budget")
        return self

    def predict(self, rows):
        return [self.coef[0] + sum(c * (row.get(k, 0.0) - self.means[k]) / self.scales[k]
                                  for c, k in zip(self.coef[1:], self.names)) for row in rows]
