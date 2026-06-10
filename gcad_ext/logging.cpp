#include <ctime>

static double
timestamp()
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec*1e-9;
}

/* Generic printing
 */

static inline void fprint(FILE *f, int n) { fprintf(f, "%d", n); }
static inline void fprint(FILE *f, long n) { fprintf(f, "%ld", n); }
static inline void fprint(FILE *f, unsigned int n) { fprintf(f, "%u", n); }
static inline void fprint(FILE *f, unsigned long n) { fprintf(f, "%lu", n); }
static inline void fprint(FILE *f, double n) { fprintf(f, "%.16e", n); }
static inline void fprint(FILE *f, bool b) { fputs(b ? "true" : "false", f); }
static inline void fprint(FILE *f, const char *s) { fputs(s, f); }

template <typename T> inline void
fprint(FILE *f, const std::vector<T> &x)
{
    size_t size = x.size();
    for (size_t i = 0; i < size; i++) {
        if (i != 0) fputs(" ", f);
        fprint(f, x[i]);
    }
}

template <typename T> struct FMT { const char *fmt; const T &val; };
template <typename T> FMT<T> fmt(const char *fmt, const T &val) { return {fmt, val}; };
template <typename T> inline void fprint(FILE *f, const FMT<T> &x) { fprintf(f, x.fmt, x.val); }

/* Logging
 */

struct Logger {
    double first_timestamp;
    double last_timestamp;
    int depth;
    bool empty_block;
};

static struct Logger LOG;

static void
LOG_start()
{
    double now = timestamp();
    LOG.first_timestamp = now;
    LOG.last_timestamp = now;
    LOG.depth = 0;
    LOG.empty_block = true;
}

static void
_next_fmt(const char *&fmt, FILE *f)
{
    for (int i = 0;; i++) {
        if (fmt[i] == '{') {
            fwrite(fmt, i, 1, f);
            fmt = fmt + i + 2;
            break;
        }
        if (fmt[i] == 0) break;
    }
}

struct _sequencehack {
    template <typename... Args> _sequencehack(Args &&...) {}
};

static void
_LOG_pr_head(double now)
{
    fprintf(stderr,
            "%.3f +%.3f ",
            now - LOG.first_timestamp,
            now - LOG.last_timestamp);
    for (int i = 0; i < LOG.depth; i++) {
        fputs("│", stderr);
    }
}

template <typename... Args> static void
_LOG_pr_body(double now, const char *fmt, const Args &...args)
{
    (void)now;
    (void)_sequencehack{(_next_fmt(fmt, stderr), fprint(stderr, args), 0)...};
    fputs(fmt, stderr);
}

static void
_LOG_pr_foot(double now)
{
    (void)now;
    fputc('\n', stderr);
    fflush(stderr);
    LOG.last_timestamp = now;
    LOG.empty_block = false;
}

template <typename... Args> static void
LOG_pr(double now, const char *fmt, const Args &...args)
{
    _LOG_pr_head(now);
    _LOG_pr_body(now, fmt, args ...);
    _LOG_pr_foot(now);
}

#define pr(...) LOG_pr(timestamp(), __VA_ARGS__)

/* Unformatted block logging
 */

static void
logline(const char *text)
{
    double now = timestamp();
    _LOG_pr_head(now);
    fputs(text, stderr);
    _LOG_pr_foot(now);
}

static double
logline_block_start(const char *text)
{
    double now = timestamp();
    _LOG_pr_head(now);
    fputs("╭", stderr);
    fputs(text, stderr);
    _LOG_pr_foot(now);
    LOG.depth++;
    LOG.empty_block = true;
    return now;
}

static void
logline_block_end(double t, const char *text)
{
    double now = timestamp();
    LOG.depth--;
    if (LOG.empty_block) fputs("\033[F\033[K", stderr);
    _LOG_pr_head(now);
    fputs(LOG.empty_block ? "-" : "╰", stderr);
    fputs(text, stderr);
    fprintf(stderr, ": %.2e", now - t);
    _LOG_pr_foot(now);
}

/* Formatted block logging
 */

struct LOG_block_info {
    char *text;
    size_t size;
    double start;
};

template <typename... Args> static void
LOG_block_start(LOG_block_info *i, const char *fmt, const Args &...args)
{
    FILE *f = open_memstream(&i->text, &i->size);
    (void)_sequencehack{(_next_fmt(fmt, f), fprint(f, args), 0)...};
    fputs(fmt, f);
    fclose(f);
    _LOG_pr_head(i->start);
    fputs("╭", stderr);
    fwrite(i->text, 1, i->size, stderr);
    _LOG_pr_foot(i->start);
    LOG.depth++;
    LOG.empty_block = true;
}

static void
LOG_block_end(const LOG_block_info *i)
{
    double now = timestamp();
    LOG.depth--;
    if (LOG.empty_block) fputs("\033[F\033[K", stderr);
    _LOG_pr_head(now);
    fputs(LOG.empty_block ? "-" : "╰", stderr);
    fwrite(i->text, 1, i->size, stderr);
    fprintf(stderr, ": %.2e", now - i->start);
    _LOG_pr_foot(now);
    free(i->text);
}

#define LOGBLOCK(...) \
    __attribute__((cleanup(LOG_block_end))) LOG_block_info _block_info = {NULL, 0, timestamp()}; \
    LOG_block_start(&_block_info, __VA_ARGS__);

#define LOGME LOGBLOCK(__func__)
