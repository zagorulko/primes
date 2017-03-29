#define _GNU_SOURCE
#include <assert.h>
#include <math.h>
#include <pthread.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>
#include <gmp.h>
#include "sieve.h"

static double elapsed_time(const struct timespec *start,
                           const struct timespec *finish)
{
    return (finish->tv_sec - start->tv_sec) +
           (finish->tv_nsec - start->tv_nsec) / 1000000000.0;
}

static void seed_rstate(gmp_randstate_t rstate, int words)
{
    FILE *urandom = fopen("/dev/urandom","r");
    assert(urandom);
    for (int i = 0; i < words; i++) {
        unsigned long seed;
        fread(&seed, sizeof(seed), 1, urandom);
        gmp_randseed_ui(rstate, seed);
    }
    fclose(urandom);
}

static inline bool try_mods(const prime_t *mods, uint64_t delta)
{
    for (int i = 1; i < SIEVE_SIZE; i++)
        if (((mods[i]+delta)%SIEVE[i]) <= 1)
            return false;
    return true;
}

static void find_probable_prime(mpz_t x, gmp_randstate_t rstate,
                                const mpz_t rand_min, const mpz_t rand_bound)
{
    prime_t mods[SIEVE_SIZE];
    uint32_t max_delta = UINT32_MAX-SIEVE[SIEVE_SIZE-1];
    bool ok = false;

    while (!ok) {
        mpz_urandomm(x, rstate, rand_bound);
        mpz_add(x, x, rand_min);
        mpz_setbit(x, 0);

        for (int i = 1; i < SIEVE_SIZE; i++)
            mods[i] = mpz_fdiv_ui(x, SIEVE[i]);

        uint32_t delta = 0;
        while (delta < max_delta) {
            if (!try_mods(mods, delta)) {
                delta += 2;
            } else {
                mpz_add_ui(x, x, delta);
                ok = true;
                break;
            }
        }
    }
}

static bool witness(mpz_t w, const mpz_t a, const mpz_t a1, const mpz_t a1_odd,
                    mp_bitcnt_t k)
{
    /* w := w^a1_odd (mod a) */
    mpz_powm(w, w, a1_odd, a);

    if (mpz_cmp_ui(w, 1) == 0)
        return false;
    /* w == -1 (mod a) */
    if (mpz_cmp(w, a1) == 0)
        return false;

    while (--k) {
        /* w := w^2 mod a */
        mpz_powm_ui(w, w, 2, a);

        if (mpz_cmp_ui(w, 1) == 0)
            return true;
        /* w == -1 (mod a) */
        if (mpz_cmp(w, a1) == 0)
            return false;
    }

    return true;
}

/*
 * Based on BN_prime_checks_for_size.
 * Should result in an error rate of less than 2^-80.
 */
static inline int checks_for_num(const mpz_t a)
{
    size_t b = mpz_sizeinbase(a, 2);
    return ((b) >= 1300 ?  2 :
            (b) >=  850 ?  3 :
            (b) >=  650 ?  4 :
            (b) >=  550 ?  5 :
            (b) >=  450 ?  6 :
            (b) >=  400 ?  7 :
            (b) >=  350 ?  8 :
            (b) >=  300 ?  9 :
            (b) >=  250 ? 12 :
            (b) >=  200 ? 15 :
            (b) >=  150 ? 18 :
            /* b >= 100 */ 27);
}

static bool miller_rabin(const mpz_t a, gmp_randstate_t rstate)
{
    assert(mpz_tstbit(a, 0));
    assert(mpz_cmp_ui(a, 1) > 0);

    /* a1 := a - 1 */
    mpz_t a1;
    mpz_init_set(a1, a);
    mpz_sub_ui(a1, a1, 1);

    /* a1 == a1_odd * 2^k */
    mp_bitcnt_t k = 1;
    while (!mpz_tstbit(a1, k))
        k++;
    mpz_t a1_odd;
    mpz_init(a1_odd);
    mpz_tdiv_q_2exp(a1_odd, a1, k);

    int checks = checks_for_num(a);
    bool ret = true;
    mpz_t w;
    mpz_init(w);

    for (int i = 0; i < checks; i++) {
        mpz_urandomm(w, rstate, a1);
        mpz_add_ui(w, w, 1);

        if (witness(w, a, a1, a1_odd, k)) {
            ret = false;
            break;
        }
    }

    mpz_clear(a1);
    mpz_clear(a1_odd);
    mpz_clear(w);
    return ret;
}

struct shared_data {
    pthread_mutex_t mutex;
    bool ok;
    mpz_t res;
    int digits;
};

struct worker_data {
    pthread_t thread;
    struct shared_data *shared;
    int i;
    double prob_time;
    double test_time;
    double mean_prob_time;
    double mean_test_time;
};

static void *worker(void *data_ptr)
{
    struct worker_data *data = data_ptr;
    struct timespec start, finish;
    double elapsed;

    gmp_randstate_t rstate;
    gmp_randinit_default(rstate);

    mpz_t rand_min, rand_bound;
    mpz_init(rand_min);
    mpz_init(rand_bound);
    mpz_ui_pow_ui(rand_min, 10, data->shared->digits-1);
    mpz_ui_pow_ui(rand_bound, 10, data->shared->digits);
    mpz_sub(rand_bound, rand_bound, rand_min);

    mpz_t x;
    mpz_init(x);

    while (!data->shared->ok) {
        if (data->i++ % 32 == 0)
            seed_rstate(rstate, 1);

        clock_gettime(CLOCK_MONOTONIC, &start);
        find_probable_prime(x, rstate, rand_min, rand_bound);
        clock_gettime(CLOCK_MONOTONIC, &finish);
        elapsed = elapsed_time(&start, &finish);
        data->prob_time += elapsed;
        data->mean_prob_time += (elapsed - data->mean_prob_time) / data->i;

        clock_gettime(CLOCK_MONOTONIC, &start);
        bool ok = miller_rabin(x, rstate);
        clock_gettime(CLOCK_MONOTONIC, &finish);
        elapsed = elapsed_time(&start, &finish);
        data->test_time += elapsed;
        data->mean_test_time += (elapsed - data->mean_test_time) / data->i;

        if (ok) {
            pthread_mutex_lock(&data->shared->mutex);
            if (!data->shared->ok) {
                mpz_set(data->shared->res, x);
                data->shared->ok = true;
            }
            pthread_mutex_unlock(&data->shared->mutex);
        }
    }

    gmp_randclear(rstate);
    mpz_clear(rand_min);
    mpz_clear(rand_bound);
    mpz_clear(x);
    return NULL;
}

void find_prime_mpz(mpz_t res, int digits)
{
    struct shared_data shared;
    pthread_mutex_init(&shared.mutex, NULL);
    shared.ok = false;
    mpz_init(shared.res);
    shared.digits = digits;

    int nthreads = sysconf(_SC_NPROCESSORS_ONLN);
    struct worker_data *threads = calloc(nthreads, sizeof(struct worker_data));
    for (int i = 0; i < nthreads; i++)
        threads[i].shared = &shared;

    struct timespec start, finish;
    clock_gettime(CLOCK_MONOTONIC, &start);
    for (int i = 0; i < nthreads; i++)
        pthread_create(&threads[i].thread, NULL, worker, &threads[i]);
    fprintf(stderr, "Searching for a %d-digit prime number using %d "
                    "threads...\n", digits, nthreads);
    for (int i = 0; i < nthreads; i++)
        pthread_join(threads[i].thread, NULL);
    clock_gettime(CLOCK_MONOTONIC, &finish);

    double total_time = elapsed_time(&start, &finish);
    int attempts = 0;
    double prob_time = 0;
    double test_time = 0;
    double mean_prob_time = 0;
    double mean_test_time = 0;
    for (int i = 0; i < nthreads; i++) {
        attempts += threads[i].i;
        prob_time += threads[i].prob_time;
        test_time += threads[i].test_time;
        mean_prob_time += threads[i].mean_prob_time;
        mean_test_time += threads[i].mean_test_time;
    }
    prob_time /= nthreads;
    test_time /= nthreads;
    mean_prob_time /= nthreads;
    mean_test_time /= nthreads;
    fprintf(stderr, "Done!\n");
    fprintf(stderr, "Total time: %f s\n", total_time);
    fprintf(stderr, "Attempts: %i\n", attempts);
    fprintf(stderr, "Probable prime search time: %f s\n", prob_time);
    fprintf(stderr, "Miller-Rabin test time: %f s\n", test_time);
    fprintf(stderr, "Mean probable prime search time: %f s\n", mean_prob_time);
    fprintf(stderr, "Mean Miller-Rabin test time: %f s\n", mean_test_time);

    // assert(mpz_probab_prime_p(shared.res, 15));
    mpz_set(res, shared.res);

    pthread_mutex_destroy(&shared.mutex);
    mpz_clear(shared.res);
    free(threads);
}

char *find_prime(int digits)
{
    mpz_t p;
    mpz_init(p);
    find_prime_mpz(p, digits);
    char *s = mpz_get_str(NULL, 10, p);
    mpz_clear(p);
    return s;
}
