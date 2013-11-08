/* Liang Xia, April 2013, 
 * modular exponentiation implementation
 * Blakley's shift-add method, binary method
 */
#include <cstdio>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <string>
#include <sys/time.h>
#include <vector>

using namespace std;

#define LEN         64
//  256 bits: LEN should be 16
//  512 bits: LEN should be 32
// 1024 bits: LEN should be 64
// 2048 bits: LEN should be 128
// 4096 bits: LEN should be 256
// 8192 bits: LEN should be 512
//16384 bits: LEN should be 1024
#define K           (LEN<<4)  // number of bits is (LEN*32 / 2)
#define MAX_UINT32  0xffffffff
#define M_ARY       8

//#define NO_SPACE
#define SHOW_ZERO

class Bignum {
    // 256 bits requires 8 elements, and 64 bits requires 2 elements.
    uint32_t num[LEN]; // little endian
public:
    static uint32_t rand_uint32(uint32_t min, uint32_t max);
    static int compare(const Bignum& b1, const Bignum& b2);
    static Bignum Blakley_shiftadd(const Bignum& a, const Bignum& b, const Bignum& n);
public:
    Bignum();
    ~Bignum() {}
    Bignum(uint32_t value);
    Bignum(const Bignum& other);
    Bignum& operator=(const Bignum& other);
    void print() const;
    int getBit(int bit) const;
    int getTotalBits() const;
    void genBignum();
    Bignum add(const Bignum& other);
    Bignum sub2(const Bignum& other); // num should be bigger than other.num
    void shiftL();
    Bignum mod_exp_binary_Blakley_shiftadd(const Bignum& exp, const Bignum& n);
};

/* start of definition of member functions */

Bignum::Bignum()
{
    memset(num, 0, sizeof(num));
}

Bignum::Bignum(uint32_t value)
{
    memset(num, 0, sizeof(num));
    num[0] = value;
}

Bignum::Bignum(const Bignum& other)
{
    memcpy(num, other.num, sizeof(num));
}

Bignum& Bignum::operator=(const Bignum& other)
{
    if (this != &other) {
        memcpy(num, other.num, sizeof(num));
    }
    return *this;
}

void Bignum::print() const
{
    for(int i = LEN - 1; i >= 0; i--) {
#ifndef SHOW_ZERO
      if (num[i] != 0) {
#endif
#ifdef NO_SPACE
        printf("%08x", num[i]);
#else
        printf("%08x ", num[i]);
#endif
#ifndef SHOW_ZERO
      }
#endif
    }
    printf("\n");
}

int Bignum::getBit(int bit) const
{
    uint32_t segment = num[bit>>5];
    bit = bit & 31;  //TODO: optimize!!
    uint32_t temp = 1 << bit;
    segment = temp & segment;
    if (temp == segment)        return 1;
    else                        return 0;
}

int Bignum::getTotalBits() const
{
    int k = LEN*32 - 1;
    while (0 == getBit(k))
        k--;
    return (k+1);
}

int Bignum::compare(const Bignum& b1, const Bignum& b2)
{
    int result = 0;
    for(int i = LEN-1; i >= 0; i--) {
        if ( b1.num[i] > b2.num[i] ) {
            result = 1;
            break;
        }
        else if (b1.num[i] < b2.num[i]) {
            result = -1;
            break;
        }
    }
    return result;
}

uint32_t Bignum::rand_uint32(uint32_t min, uint32_t max)
{
    uint32_t result = 0;
    if (min > max) {
        printf("Bignum::gen wrong min and max\n");
        return result;
    }
    uint64_t range = static_cast<uint64_t>(max) - static_cast<uint64_t>(min) + 1;
    result = (uint32_t) ( ( rand() % range ) + min ); //TODO
    return result;
}

void Bignum::genBignum()
{
    for (int i = 0; i < LEN >> 1; i++) {
        num[i] = rand_uint32(0, MAX_UINT32);
    }
}

Bignum Bignum::add(const Bignum& other)
{
    Bignum result;
    uint64_t temp;
    uint64_t carry = 0;

    for(int i = 0; i < LEN; i++) {
        temp = static_cast<uint64_t>(num[i]) 
               + static_cast<uint64_t>(other.num[i])
               + carry;
        carry = temp >> 32;
        result.num[i] = temp & MAX_UINT32;
    }
    return result;
}

// sub2 will happen only if this->num is bigger than other.num
Bignum Bignum::sub2(const Bignum& other)
{
    Bignum result;
    uint64_t temp;
    uint64_t carry = 0;

    for (int i = 0; i < LEN; i++) {
        uint64_t num_a = static_cast<uint64_t>(this->num[i]);
        uint64_t num_b = static_cast<uint64_t>(other.num[i]);
        if ( num_a >= num_b + carry) {
            temp = num_a - num_b - carry;
            carry = 0;
        }
        else {
            temp = num_a + MAX_UINT32 + 1 - num_b - carry;
            carry = 1;
        }
        result.num[i] = static_cast<uint32_t>(temp);
    }
    return result;
}

void Bignum::shiftL()
{
    uint32_t CONSTANT = 0x80000000;
    uint32_t carry = 0;
    for (int i = 0; i < LEN; i++) {
        uint32_t temp = num[i];
        num[i] = (num[i] << 1) + carry;
        temp = temp & CONSTANT;
        if (temp == CONSTANT) {
            carry = 1;
        }
        else {
            carry = 0;
        }
    }
}

Bignum Bignum::mod_exp_binary_Blakley_shiftadd(const Bignum& exp, const Bignum& n)
{
    Bignum C(1);
    Bignum M = *this;
    int k = exp.getTotalBits();
    if (1 == exp.getBit(k-1))
        C = M;
    for (int i = k-2; i >= 0; i--) {
        C = Blakley_shiftadd(C, C, n);
        if (1 == exp.getBit(i)) {
            C = Blakley_shiftadd(C, M, n);
        }
    }
    return C;
}

Bignum Bignum::Blakley_shiftadd(const Bignum& a, const Bignum& b, const Bignum& n)
{
    Bignum R;
    for (int i = 0; i < K; i++) {
        R.shiftL();
        if ( a.getBit(K-1-i) == 1)
            R = R.add(b);
        while (compare(R, n) >= 0)
            R = R.sub2(n);
    }
    return R;
}

/* end of definition of member functions */

/* start of definition of local functions */

double read_timer()
{
    static bool initialized = false;
    static struct timeval start;
    struct timeval end;
    if( !initialized )
    {
        gettimeofday( &start, NULL );
        initialized = true;
    }

    gettimeofday( &end, NULL );

    return (end.tv_sec - start.tv_sec) + 1.0e-6 * (end.tv_usec - start.tv_usec);
}

void test8()
{
    srand (time(NULL));

    Bignum M;
    M.genBignum();
    printf("  M = "); M.print(); printf("\n");

    Bignum exp;
    exp.genBignum();
    printf("exp = "); exp.print(); printf("\n");

    Bignum n;
    n.genBignum();
    printf("  n = "); n.print(); printf("\n");

    double seconds;
    Bignum m1, m2;

    printf("Multiplication optimization. \n");
    printf("           exponentiation - binary method\n");
    printf("           multiplication - shift-add (Blakley method)\n\n");
    seconds = read_timer();
    Bignum re3 = M.mod_exp_binary_Blakley_shiftadd(exp, n);
    seconds = read_timer() - seconds;
    printf("re3 = "); re3.print(); printf("\n");
    printf("\ntime = %lf\n\n", seconds);
}


int main(int argc, char** argv)
{
    printf("bit sizes = %d bits\n\n", K);
    test8();
    return 0;
}

