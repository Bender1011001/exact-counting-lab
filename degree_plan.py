"""Certified occurrence-wise degree bound and sufficient Hermite conditions."""


def form_degree_bound(design):
    """deg(prod_e Phi_e(P_e(T))) <= sum_(e,k) deg(P_e,k).

    Zero and constant forms contribute zero. This bound ignores cancellation
    between factors and therefore remains safe without knowing the count.
    """
    return sum(max((j for coefficients in form.values()
                    for j, c in enumerate(coefficients) if c != 0), default=0)
               for clause in design.P for form in clause)


def recovery_plan(design, max_order):
    if design.Q < 1 or max_order < 0:
        raise ValueError('Positive point count and nonnegative order required')
    degree = form_degree_bound(design)
    base, extra = divmod(degree+1, design.Q)
    orders = {t: base + int(t <= extra) - 1 for t in range(1,design.Q+1)
              if base + int(t <= extra)}
    if max(orders.values()) > max_order:
        raise ValueError('Requested maximum Taylor order cannot certify recovery')
    return degree, orders
