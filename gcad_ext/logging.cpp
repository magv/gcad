/* ╷╷     ╷
 * ├┤╭╮╭┐ │ hierarchical
 * ╵╵╰╯╰┤ │ log
 *     ╶╯ ╵
 */

#include <atomic>
#include <mutex>
#include <stdio.h>
#include <thread>
#include <time.h>
#include <unistd.h>
#include <vector>

/* Generic formatting
 */

static inline void fprint(FILE *f, int n) { fprintf(f, "%d", n); }
static inline void fprint(FILE *f, long n) { fprintf(f, "%ld", n); }
static inline void fprint(FILE *f, long long n) { fprintf(f, "%lld", n); }
static inline void fprint(FILE *f, unsigned int n) { fprintf(f, "%u", n); }
static inline void fprint(FILE *f, unsigned long n) { fprintf(f, "%lu", n); }
static inline void fprint(FILE *f, unsigned long long n) { fprintf(f, "%llu", n); }
static inline void fprint(FILE *f, float n) { fprintf(f, "%.8e", n); }
static inline void fprint(FILE *f, double n) { fprintf(f, "%.16e", n); }
static inline void fprint(FILE *f, bool b) { fputs(b ? "true" : "false", f); }
static inline void fprint(FILE *f, const char *s) { fputs(s, f); }

template <typename T> static void
fprint(FILE *f, const std::vector<T> &x)
{
    size_t size = x.size();
    for (size_t i = 0; i < size; i++) {
        if (i != 0) fputs(" ", f);
        fprint(f, x[i]);
    }
}

template <typename T> struct FMT { const char *fmt; const T &val; };
template <typename T> static FMT<T> fmt(const char *fmt, const T &val) { return {fmt, val}; };
template <typename T> static inline void fprint(FILE *f, const FMT<T> &x) { fprintf(f, x.fmt, x.val); }

struct FMT_SEC { const double val; };
static inline FMT_SEC fmt_sec(double sec) { return {sec}; };
static inline void fprint(FILE *f, const FMT_SEC x) {
    const char *unit;
    double scale;
    if (x.val < 1e-6) { unit = "ns"; scale = 1e9; }
    else if (x.val < 1e-3) { unit = "µs"; scale = 1e6; }
    else if (x.val < 1.) { unit = "ms"; scale = 1e3; }
    else if (x.val < 60.) { unit = "s"; scale = 1.; }
    else if (x.val < 3600.) { unit = "m"; scale = 1./60.; }
    else if (x.val < 24.*3600.) { unit = "h"; scale = 1./3600.; }
    else if (x.val < 7.*24.*3600.) { unit = "d"; scale = 1./(24.*3600.); }
    else if (x.val < 365.*24.*3600.) { unit = "w"; scale = 1./(7.*24.*3600.); }
    else { unit = "y"; scale = 1./(365.*24.*3600.); };
    fprintf(f, "%.3f%s", x.val*scale, unit);
}

struct FMT_B { const double val; };
static inline FMT_B fmt_bytes(double sec) { return {sec}; };
static inline void fprint(FILE *f, const FMT_B x) {
    const char *unit;
    double scale;
    if (x.val < (1<<10)) { unit = "B"; scale = 1.; }
    else if (x.val < (1ul<<20)) { unit = "kB"; scale = 1./(1ul<<10); }
    else if (x.val < (1ul<<30)) { unit = "MB"; scale = 1./(1ul<<20); }
    else if (x.val < (1ul<<40)) { unit = "GB"; scale = 1./(1ul<<30); }
    else if (x.val < (1ul<<50)) { unit = "TB"; scale = 1./(1ul<<40); }
    else { unit = "PB"; scale = 1./(1ul<<50); }
    fprintf(f, "%.1f%s", x.val*scale, unit);
}

template <typename... Args> static void
fprint_fmt(FILE *f, const char *fmt, const Args &...args)
{
    auto next_format = [](const char *fmt, ssize_t idx) {
        for (ssize_t i = idx;; i++) {
            if (fmt[i] == '{') return i + 2;
            if (fmt[i] == 0) return i;
        }
    };
    (void)next_format;
    ssize_t idx = 0, next;
    ((next = next_format(fmt, idx),
      (next - idx - 2 > 0) ? fwrite(fmt + idx, next - idx - 2, 1, f) : 0,
      idx = next,
      fprint(f, args)),
     ...);
    fputs(fmt + idx, f);
}

/* Logging base
 */

using span_t = int;

struct Span {
    span_t parent = -1;
    std::vector<span_t> children;
    char *text = NULL;
    size_t text_size = 0;
    double ndone = 0;
    double ntotal = 0;
    bool active = false;
    double start_time = 0;
    int depth = 0;
    int log_depth = 0;
};

static double
timestamp()
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec * 1e-9;
}

static void
_print_progress_bar(FILE *out, double frac, int width)
{
    static const char *fraction[] = {"╶", "╸"};
    double f = width * frac;
    int n = (int)f;
    if (f < 0.5) {
        for (int i = 0; i < width; i++)
            fputs("─", out);
    } else {
        for (int i = 0; i < n; i++)
            fputs("━", out);
        if (n < width) {
            fputs(fraction[int((f - n) * 2)], out);
            for (int i = n + 1; i < width; i++)
                fputs("─", out);
        }
    }
}

static void
_print_tree_prefix(FILE *out, int depth, const int *is_last_ancestors)
{
    if (depth == 0) {
        fputs("⯁ ", out);
    } else {
        fputs("  ", out);
        for (int i = 1; i < depth; i++) {
            fputs(is_last_ancestors[i] ? "   " : "│  ", out);
        }
        fputs(is_last_ancestors[depth] ? "╰─ " : "├─ ", out);
    }
}

struct Hog {
    std::vector<Span> spans;
    std::vector<span_t> free_list;
    std::vector<int> is_last_ancestors;
    int last_line_count = 0;
    FILE *out;
    std::atomic<bool> is_running{false};
    std::thread render_thread;
    std::mutex mtx;
    std::atomic<bool> must_repaint{false};

    Hog(FILE *f = stderr) : out(f)
    {
        Span root;
        root.depth = -1;
        root.log_depth = -1;
        root.active = true;
        root.text = NULL;
        root.text_size = 0;
        root.start_time = timestamp();
        spans.push_back(root);
    }

    ~Hog()
    {
        is_running.store(false);
        if (render_thread.joinable()) render_thread.join();
    }

    void
    start()
    {
        is_running.store(true);
        render_thread = std::thread([this]() {
            while (is_running.load()) {
                usleep(200000);
                if (must_repaint.exchange(false)) {
                    std::lock_guard<std::mutex> lock(mtx);
                    render();
                }
            }
        });
    }

    void
    stop()
    {
        is_running.store(false);
        must_repaint.store(false);
        if (render_thread.joinable()) render_thread.join();
        fprintf(out, "\033[J");
        fflush(out);
        last_line_count = 0;
    }

    template <typename... Args> void
    log(span_t parent, double now, const char *fmt, const Args &...args)
    {
        std::lock_guard<std::mutex> lock(mtx);
        double elapsed = now - spans[0].start_time;
        fprintf(out, "\033[2m%5.3f\033[0m ", elapsed);
        int d = spans[parent].log_depth + 1;
        for (int i = 0; i < d; i++) {
            fputs("│ ", out);
        }
        fprint_fmt(out, fmt, args...);
        fputs("\033[K\n\033[K", out);
        must_repaint.store(false);
        render();
    }

    span_t
    _alloc_span(span_t parent, double now, char *text, size_t size)
    {
        Span s;
        s.parent = parent;
        s.text = text;
        s.text_size = size;
        s.depth = spans[parent].depth + 1;
        s.log_depth = spans[parent].log_depth;
        s.active = true;
        s.start_time = now;
        span_t id;
        if (free_list.empty()) {
            id = (span_t)spans.size();
            spans.push_back(s);
        } else {
            id = free_list.back();
            free_list.pop_back();
            spans[id] = s;
        }
        spans[parent].children.push_back(id);
        return id;
    }

    template <typename... Args> span_t
    begin_span(span_t parent,
               double now,
               bool print_log,
               const char *fmt,
               const Args &...args)
    {
        char *text = NULL;
        size_t size = 0;
        FILE *f = open_memstream(&text, &size);
        fprint_fmt(f, fmt, args...);
        fclose(f);
        span_t id = _alloc_span(parent, now, text, size);
        if (print_log) {
            spans[id].log_depth += 1;
            log(parent, now, "╭ {}", text);
        } else {
            must_repaint.store(true);
        }
        return id;
    }

    span_t
    end_span(span_t id, double now, bool print_log)
    {
        auto &s = spans[id];
        span_t parent = s.parent;
        double start = s.start_time;
        s.active = false;
        auto &children = spans[parent].children;
        for (auto it = children.begin(); it != children.end(); ++it) {
            if (*it == id) {
                children.erase(it);
                break;
            }
        }
        s.parent = -1;
        s.children.clear();
        free_list.push_back(id);
        if (print_log) {
            log(parent, now, "╰ {} \033[2m({})\033[0m", s.text, fmt_sec(now - start));
        } else {
            must_repaint.store(true);
        }
        return parent;
    }

    template <typename... Args> void
    set_text(span_t id, const char *fmt, const Args &...args)
    {
        Span &s = spans[id];
        FILE *f = open_memstream(&s.text, &s.text_size);
        fprint_fmt(f, fmt, args ...);
        fclose(f);
        must_repaint.store(true);
    }

    void
    set_progress(span_t id, double done, double total)
    {
        spans[id].ndone = done;
        spans[id].ntotal = total;
        must_repaint.store(true);
    }

    int
    _render_span(span_t id, int depth, int *is_last_ancestors)
    {
        auto &s = spans[id];
        _print_tree_prefix(out, depth, is_last_ancestors);
        fputs(s.text, out);
        if (s.ntotal > 0) {
            double frac = s.ndone / s.ntotal;
            fputs(" ", out);
            _print_progress_bar(out, frac, 20);
        }
        fputs("\033[K\n", out);
        int count = 1;
        if (!s.children.empty()) {
            int nchildren = (int)s.children.size();
            for (int vi = 0; vi < nchildren; vi++) {
                is_last_ancestors[depth + 1] = vi + 1 == nchildren;
                count += _render_span(s.children[vi], depth + 1, is_last_ancestors);
            }
        }
        return count;
    }

    void
    render()
    {
        if (is_last_ancestors.size() < spans.size()) {
            is_last_ancestors.resize(spans.size());
        }
        int line_count = 1;
        fputs("\n", out);
        int nchildren = (int)spans[0].children.size();
        for (int i = 0; i < nchildren; i++) {
            is_last_ancestors[0] = i + 1 == nchildren;
            line_count += _render_span(spans[0].children[i], 0, is_last_ancestors.data());
        }
        fputs("\033[J", out);
        for (; line_count < last_line_count; line_count++) {
            fputs("\n", out);
        }
        fprintf(out, "\033[%dA", line_count);
        fflush(out);
        last_line_count = line_count;
    }
};

/* Logging & tracing
 */

static Hog HOG;
static span_t HOG_current_span = 0;

#define log_start() HOG.start()

#define log_stop() HOG.stop()

#define log(...) HOG.log(HOG_current_span, timestamp(), __VA_ARGS__)

#define trace_enter(...) \
    do { HOG_current_span = HOG.begin_span(HOG_current_span, timestamp(), false, __VA_ARGS__); } while(0)

#define trace_exit() \
    do { HOG_current_span = HOG.end_span(HOG_current_span, timestamp(), false); } while(0)

#define log_trace_enter(...) \
    do { HOG_current_span = HOG.begin_span(HOG_current_span, timestamp(), true, __VA_ARGS__); } while(0)

#define log_trace_exit() \
    do { HOG_current_span = HOG.end_span(HOG_current_span, timestamp(), true); } while(0)

#define trace_text(...) HOG.set_text(HOG_current_span, __VA_ARGS__)

#define trace_progress(ndone, ntotal) HOG.set_progress(HOG_current_span, ndone, ntotal)

/* Scope-based logging & tracing
 */

struct _HOG_logged_scope_info {};
struct _HOG_silent_scope_info {};

static void
_HOG_logged_scope_exit(const _HOG_logged_scope_info *i)
{
    (void)i;
    log_trace_exit();
}

static void
_HOG_silent_scope_exit(const _HOG_silent_scope_info *i)
{
    (void)i;
    trace_exit();
}

#define log_trace_scope(...) \
    __attribute__((cleanup(_HOG_logged_scope_exit))) _HOG_logged_scope_info _scope_info = {}; \
    log_trace_enter(__VA_ARGS__);

#define trace_scope(...) \
    __attribute__((cleanup(_HOG_silent_scope_exit))) _HOG_silent_scope_info _scope_info = {}; \
    trace_enter(__VA_ARGS__);
