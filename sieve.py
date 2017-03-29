#!/usr/bin/env python3
import argparse
import math

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('n', type=int)
    args = parser.parse_args()

    primes = [2]
    p = 1
    while len(primes) < args.n:
        p += 2
        s = int(math.sqrt(p))
        ok = True
        for q in primes:
            if q > s:
                break
            if p % q == 0:
                ok = False
                break
        if ok:
            primes.append(p)

    with open('sieve.h', 'w') as fout:
        fout.write('#define SIEVE_SIZE {}\n'.format(args.n))
        fout.write('typedef uint32_t prime_t;\n')
        fout.write('static const prime_t SIEVE[SIEVE_SIZE] = {\n')
        for p in primes:
            fout.write('{},\n'.format(p))
        fout.write('};\n')

if __name__ == '__main__':
    main()
