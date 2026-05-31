/* The following is a minimal wrapper around FLINT's types that
 * will help our code to represent the underlying math tersely,
 * with no loss of performance.
 */

#include <flint/flint.h>
#include <flint/fmpq.h>
#include <flint/fmpz.h>
#include <flint/fmpz_poly.h>
#include <flint/fmpz_mpoly.h>

struct Q;
struct Z;

struct Z {
    fmpz z;
    Z() { fmpz_init(&z); }
    Z(Z &&x) { z = x.z; fmpz_init(&x.z); }
    Z(const Z &x) { fmpz_init_set(&z, &x.z); }
    Z &operator =(Z &&x) { fmpz_swap(&z, &x.z); return *this; }
    Z &operator =(const Z &x) { fmpz_set(&z, &x.z); return *this; }
    ~Z() { fmpz_clear(&z); }
    bool is_zero() const { return fmpz_is_zero(&z); }
    bool is_one() const { return fmpz_is_one(&z); }
    int sign() const { return fmpz_sgn(&z); }
    bool is_positive() const { return sign() == +1; }
    bool is_negative() const { return sign() == -1; }
    double get_d() const { return fmpz_get_d(&z); }
    Z &operator *=(const Z &x) { fmpz_mul(&z, &z, &x.z); return *this; }
    Z &operator +=(const Z &x) { fmpz_add(&z, &z, &x.z); return *this; }
    Z &operator -=(const Z &x) { fmpz_sub(&z, &z, &x.z); return *this; }
    bool operator <(const Z &x) const { return fmpz_cmp(&z, &x.z) < 0; }
    bool operator <=(const Z &x) const { return fmpz_cmp(&z, &x.z) <= 0; }
    bool operator ==(const Z &x) const { return fmpz_cmp(&z, &x.z) == 0; }
    bool operator >=(const Z &x) const { return fmpz_cmp(&z, &x.z) >= 0; }
    bool operator >(const Z &x) const { return fmpz_cmp(&z, &x.z) > 0; }
    void set_neg() { fmpz_neg(&z, &z); }
    Z pow_ui(ulong e) const { Z r; fmpz_pow_ui(&r.z, &z, e); return r; }
};

struct Q {
    fmpq q;
    Q() { fmpq_init(&q); }
    Q(Q &&x) { q = x.q; fmpq_init(&x.q); }
    Q(const Q &x) { fmpq_init(&q); fmpq_set(&q, &x.q); }
    Q &operator =(Q &&x) { fmpq_swap(&q, &x.q); return *this; }
    Q &operator =(const Q &x) { fmpq_set(&q, &x.q); return *this; }
    ~Q() { fmpq_clear(&q); }
    bool is_zero() const { return fmpq_is_zero(&q); }
    bool is_one() const { return fmpq_is_one(&q); }
    int sign() const { return fmpq_sgn(&q); }
    bool is_positive() const { return sign() == +1; }
    bool is_negative() const { return sign() == -1; }
    double get_d() const { return fmpq_get_d(&q); }
    Q &operator *=(const Q &x) { fmpq_mul(&q, &q, &x.q); return *this; }
    Q &operator +=(const Q &x) { fmpq_add(&q, &q, &x.q); return *this; }
    Q &operator -=(const Q &x) { fmpq_sub(&q, &q, &x.q); return *this; }
    Q &operator /=(const Q &x) { fmpq_div(&q, &q, &x.q); return *this; }
    bool operator <(const Q &x) const { return fmpq_cmp(&q, &x.q) < 0; }
    bool operator <=(const Q &x) const { return fmpq_cmp(&q, &x.q) <= 0; }
    bool operator ==(const Q &x) const { return fmpq_cmp(&q, &x.q) == 0; }
    bool operator !=(const Q &x) const { return fmpq_cmp(&q, &x.q) != 0; }
    bool operator >=(const Q &x) const { return fmpq_cmp(&q, &x.q) >= 0; }
    bool operator >(const Q &x) const { return fmpq_cmp(&q, &x.q) > 0; }
    void set_inv() { fmpz_swap(fmpq_numref(&q), fmpq_denref(&q)); }
    void set_mul_2exp(int e) { if (e >= 0) { fmpq_mul_2exp(&q, &q, e); } else { fmpq_div_2exp(&q, &q, -e); } }
    void set_neg() { fmpq_neg(&q, &q); }
    Q pow_si(slong e) const { Q r; fmpq_pow_si(&r.q, &q, e); return r; }
    Q inv() const { Q r; fmpz_set(fmpq_numref(&r.q), fmpq_denref(&q)); fmpz_set(fmpq_denref(&r.q), fmpq_numref(&q)); return r; }
    Z truncate() const { Z r; fmpz_tdiv_q(&r.z, fmpq_numref(&q), fmpq_denref(&q)); return r; }
};

static Z operator +(const Z &a, const Z &b) { Z r; fmpz_add(&r.z, &a.z, &b.z); return r; }
static Z operator -(const Z &a, const Z &b) { Z r; fmpz_sub(&r.z, &a.z, &b.z); return r; }
static Z operator *(const Z &a, const Z &b) { Z r; fmpz_mul(&r.z, &a.z, &b.z); return r; }
static Q operator /(const Z &a, const Z &b) { Q r; fmpq_set_fmpz_frac(&r.q, &a.z, &b.z); return r; }
static Z operator -(const Z &a) { Z r; fmpz_neg(&r.z, &a.z); return r; }

static Q operator +(const Q &a, const Q &b) { Q r; fmpq_add(&r.q, &a.q, &b.q); return r; }
static Q operator -(const Q &a, const Q &b) { Q r; fmpq_sub(&r.q, &a.q, &b.q); return r; }
static Q operator *(const Q &a, const Q &b) { Q r; fmpq_mul(&r.q, &a.q, &b.q); return r; }
static Q operator /(const Q &a, const Q &b) { Q r; fmpq_div(&r.q, &a.q, &b.q); return r; }
static Q operator -(const Q &a) { Q r; fmpq_neg(&r.q, &a.q); return r; }

struct ZPoly {
    fmpz_poly_struct p;
    ZPoly() { fmpz_poly_init(&p); }
    ZPoly(ZPoly &&x) { p = x.p; fmpz_poly_init(&x.p); }
    ZPoly(const ZPoly &x) { fmpz_poly_init(&p); fmpz_poly_set(&p, &x.p); }
    ZPoly &operator =(ZPoly &&x) { fmpz_poly_swap(&p, &x.p); return *this; }
    ZPoly &operator =(const ZPoly &x) { fmpz_poly_set(&p, &x.p); return *this; }
    ~ZPoly() { fmpz_poly_clear(&p); }
    slong length() const { return fmpz_poly_length(&p); }
    Z &operator[](size_t i) { return *(Z*)fmpz_poly_get_coeff_ptr(&p, i); }
    const Z &operator[](size_t i) const { return *(Z*)fmpz_poly_get_coeff_ptr(&p, i); }
    Z &lead() { return *(Z*)fmpz_poly_lead(&p); }
    const Z &lead() const { return *(Z*)fmpz_poly_lead(&p); }
    void set_taylor_shift(const Z &c) { fmpz_poly_taylor_shift(&p, &p, &c.z); }
    void set_reverse(slong n) { fmpz_poly_reverse(&p, &p, n); }
    void set_neg() { fmpz_poly_neg(&p, &p); }
    void set_shift_right(slong n) { fmpz_poly_shift_right(&p, &p, n); }
    void set_divexact_root(const Q &c) { fmpz_poly_divexact_root_fmpq(&p, &p, &c.q); }
    Q eval(const Q &a) const { Q r; fmpz_poly_evaluate_fmpq(&r.q, &p, &a.q); return r; }
};

#define ZPoly_of_fmpz_poly(x) (*(ZPoly*)(x))
#define ZPoly_to_fmpz_poly(x) (&(x).p)
#define Z_of_fmpz(x) (*(Z*)(x))
#define Z_to_fmpz(x) (&(x).z)
#define Q_of_fmpq(x) (*(Q*)(x))
#define Q_to_fmpq(x) (&(x).q)

static Z Z_of_si(slong x) { Z z; fmpz_set_si(&z.z, x); return z; }
static Q Q_of_si(slong n, slong d = 1) { Q q; fmpq_set_si(&q.q, n, d); return q; }
static Q Q_ldexp_si(slong m, int e) { Q q = Q_of_si(m); if (e >= 0) {fmpq_mul_2exp(&q.q, &q.q, e);} else {fmpq_div_2exp(&q.q, &q.q, -e);} return q; }
