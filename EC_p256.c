/*EC: y^2 = x^3 + A * x + B,
implemented in 2^30 base
every number is taken using 10 limbs
p[10] = {1073741823, 1073741823, 1073741823, 63, 0, 0, 4096, 1073725440, 65535, 0}
is the representation of the prime in 2^30 base*/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

void add(long long int *, long long int *, long long int *);//addition of two 10 limb numbers
void add_fp(long long int *, long long int *, long long int *);//modular addition
void subt(long long int *, long long int *, long long int *, int);//subtraction of two variable limb numbers. can only subtract smaller number from bigger number
void mult(long long int *, long long int *, long long int *);//karatsuba multiplication of two 10 limb numbers
void mult_fp(long long int *, long long int *, long long int *);//modular multiplication
void karatsuba_2limb(long long int *, long long int *, long long int *);//karatsuba multiplication for 2 limbs
void karatsuba_4limb(long long int *, long long int *, long long int *);//karatsuba multiplication for 4 limbs
void karatsuba_8limb(long long int *, long long int *, long long int *);//karatsuba multiplication of 8 limbs
void barrett(long long int *, long long int *);// barrett reduction
void inverse(long long int *, long long int *);
void pt_add(long long int *, long long int *, long long int *, long long int *, long long int *, long long int *);//elliptic curve addition of two points
void pt_double(long long int *, long long int *, long long int *, long long int *);//doubling of a point on elliptic curve
void scalar_mult(long long int *, long long int *, long long int *, long long int *, long long int *);//scalar multiplication of points on elliptic curve
int compare(long long int *, long long int *);//compare two 10 limb numbers
int first_1(long long int *, int *);//first one in the binary representation of a 9 limb number

void add(long long int *c, long long int *d, long long int *e)
{
    long long int c1[9], d1[9];
    int i;
    for (i = 0; i < 9; i++)
    {
        c1[i] = c[i];//copying the inputs so that it doesn't get changed during operation
        d1[i] = d[i];
    }
    for (int i = 0; i < 9; i++)
    {
        e[i] += (c1[i] + d1[i]);//adding corresponding limbs
        e[i + 1] = (e[i] >> 30) & 1;//checking for carry
        e[i] = e[i] & 0x3fffffff;//putting last 30 bits in the result
    }
}

void add_fp(long long int *c, long long int *d, long long int *e)
{
    long long int e1[18] = {0};
    add(c, d, e1);
    barrett(e1, e);
}

void subt(long long int *c, long long int *d, long long int *e, int limbs)//number of limbs to represent the bigger input is taken on "limbs"
{
    long long int c1[limbs], d1[limbs];
    int i;
    for (i = 0; i < limbs; i++)
    {
        c1[i] = c[i];//copying the inputs so that it doesn't get changed during operation
        d1[i] = d[i];
    }
    for (i = 0; i < limbs; i++)
    {
        e[i] = c1[i] - d1[i];//subtracting corresponding limbs
        if ((e[i] >> 63) & 1)//checking if the corresponding limb negative
        {
            e[i] += (1 << 30);//if corresponding limb is negative adding 2^30 with it
            c1[i + 1] -= 1;//and taking 1 borrow from next limb
        }
    }
}

void mult(long long int *a, long long int *b, long long int *c)
{
    int i;
    long long int a0[2] = {0}, b0[2] = {0}, c0[16] = {0};
    long long int a1_dash[8] = {0}, b1_dash[8] = {0}, c1_dash[16] = {0};
    long long int a2[8] = {0}, b2[8] = {0}, c2[16] = {0};
    long long int c1_dash_dash[16] = {0}, c1[16] = {0};
    for (i = 0; i < 2; i++)
    {
        a0[i] = a[i];
        b0[i] = b[i];
    }
    for (i = 2; i < 10; i++)
    {
        a2[i - 2] = a[i];
        b2[i - 2] = b[i];
    }
    for (i = 0; i < 2; i++)
    {
        a1_dash[i] = a0[i] + a2[i];
        b1_dash[i] = b0[i] + b2[i];
    }
    for (i = 2; i < 8; i++)
    {
        a1_dash[i] = a2[i];
        b1_dash[i] = b2[i];
    }
    karatsuba_2limb(a0, b0, c0);
    karatsuba_8limb(a1_dash, b1_dash, c1_dash);
    karatsuba_8limb(a2, b2, c2);
    subt(c1_dash, c0, c1_dash_dash, 16);
    subt(c1_dash_dash, c2, c1, 16);
    c[0] = c0[0];
    c[1] = c0[1];
    c[2] = c0[2] + c1[0];
    c[3] = c0[3] + c1[1];
    c[4] = c1[2] + c2[0];
    c[5] = c1[3] + c2[1];
    c[6] = c1[4] + c2[2];
    c[7] = c1[5] + c2[3];
    c[8] = c1[6] + c2[4];
    c[9] = c1[7] + c2[5];
    c[10] = c1[8] + c2[6];
    c[11] = c1[9] + c2[7];
    c[12] = c1[10] + c2[8];
    c[13] = c1[11] + c2[9];
    c[14] = c1[12] + c2[10];
    c[15] = c1[13] + c2[11];
    c[16] = c1[14] + c2[12];
    c[17] = c1[15] + c2[13];
    c[18] = c2[14];
    c[19] = c2[15];
    for (i = 0; i < 20; i++)
    {
        c[i + 1] += (c[i] >> 30);
        c[i] = c[i] & 0x3fffffff;
    }
}

void mult_fp(long long int *c, long long int *d, long long int *e)
{
    int i, j, k;
    long long int e1[20] = {0};
    mult(c, d, e1);
    barrett(e1, e);
}

void karatsuba_2limb(long long int *a, long long int *b, long long int *c)
{
    int i;
    c[0] = a[0] * b[0];
    c[2] = a[1] * b[1];
    c[1] = ((a[0] + a[1]) * (b[0] + b[1])) - c[0] - c[2];
}

void karatsuba_4limb(long long int *a, long long int *b, long long int *c)
{
    long long int a0[2] = {0}, b0[2] = {0}, c0[4] = {0};
    long long int a1_dash[2] = {0}, b1_dash[2] = {0}, c1_dash[4] = {0};
    long long int a2[2] = {0}, b2[2] = {0}, c2[4] = {0};
    long long int c1_dash_dash[4] = {0}, c1[4] = {0};
    int i;
    for (i = 0; i < 2; i++)
    {
        a0[i] = a[i];
        b0[i] = b[i];
    }
    for (i = 2; i < 4; i++)
    {
        a2[i - 2] = a[i];
        b2[i - 2] = b[i];
    }
    for (i = 0; i < 2; i++)
    {
        a1_dash[i] = a0[i] + a2[i];
        b1_dash[i] = b0[i] + b2[i];
    }
    karatsuba_2limb(a0, b0, c0);
    karatsuba_2limb(a1_dash, b1_dash, c1_dash);
    karatsuba_2limb(a2, b2, c2);
    subt(c1_dash, c0, c1_dash_dash, 4);
    subt(c1_dash_dash, c2, c1, 4);
    c[0] = c0[0];
    c[1] = c0[1];
    c[2] = c0[2] + c1[0];
    c[3] = c0[3] + c1[1];
    c[4] = c1[2] + c2[0];
    c[5] = c1[3] + c2[1];
    c[6] = c2[2];
    c[7] = c2[3];
}

void karatsuba_8limb(long long int *a, long long int *b, long long int *c)
{
    int i;
    long long int a0[4] = {0}, b0[4] = {0}, c0[8] = {0};
    long long int a1_dash[4] = {0}, b1_dash[4] = {0}, c1_dash[8] = {0};
    long long int a2[4] = {0}, b2[4] = {0}, c2[8] = {0};
    long long int c1_dash_dash[8] = {0}, c1[8] = {0};
    for (i = 0; i < 4; i++)
    {
        a0[i] = a[i];
        b0[i] = b[i];
    }
    for (i = 4; i < 8; i++)
    {
        a2[i - 4] = a[i];
        b2[i - 4] = b[i];
    }
    for (i = 0; i < 4; i++)
    {
        a1_dash[i] = a0[i] + a2[i];
        b1_dash[i] = b0[i] + b2[i];
    }
    karatsuba_4limb(a0, b0, c0);
    karatsuba_4limb(a1_dash, b1_dash, c1_dash);
    karatsuba_4limb(a2, b2, c2);
    subt(c1_dash, c0, c1_dash_dash, 8);
    subt(c1_dash_dash, c2, c1, 8);
    c[0] = c0[0];
    c[1] = c0[1];
    c[2] = c0[2];
    c[3] = c0[3];
    c[4] = c0[4] + c1[0];
    c[5] = c0[5] + c1[1];
    c[6] = c0[6] + c1[2];
    c[7] = c0[7] + c1[3];
    c[8] = c1[4] + c2[0];
    c[9] = c1[5] + c2[1];
    c[10] = c1[6] + c2[2];
    c[11] = c1[7] + c2[3];
    c[12] = c2[4];
    c[13] = c2[5];
    c[14] = c2[6];
    c[15] = c2[7];
}

void barrett(long long int *x, long long int *reduced_x)
{
    for (int i = 0; i < 10; i++)
        reduced_x[i] = 0;
    int i;
    long long int p[10] = {1073741823, 1073741823, 1073741823, 63, 0, 0, 4096, 1073725440, 65535, 0};
    long long int T[10] = {805306368, 0, 0, 1073741820, 1073741807, 1073741759, 1073741567, 1073741823, 4095, 16384};//precomputed T value
    long long int q0[10] = {0}, q1[20] = {0}, q2[10] = {0}, qp[20] = {0}, r[11] = {0}, r1[11] = {0}, r2[11] = {0};
    long long int arr[11] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1};
    for (i = 0; i < 10; i++)
        q0[i] = x[i + 8];
    mult(q0, T, q1);
    for (i = 0; i < 10; i++)
        q2[i] = q1[i + 10];
    mult(q2, p, qp);
    if (compare(x, qp) == 1)
        subt(x, qp, r1, 10);
    else
    {
        subt(qp, x, r, 10);
        subt(arr, r, r1, 11);
    }
    while (compare(r1, p) == 1)
    {
        subt(r1, p, r2, 10);
        for (int j = 0; j < 10; j++)
        {
            r1[j] = r2[j];
            r2[j] = 0;
        }
    }
    for (int j = 0; j < 10; j++)
        reduced_x[j] = r1[j];
}

void inverse(long long int *a, long long int *a_inv)
{
    int i, j, k;
    long long int p_minus_2[9] = {1073741821, 1073741823, 1073741823, 63, 0, 0, 4096, 1073725440, 65535};//(p-2)
    long long int x1[10] = {0}, x[10] = {1, 0};
    for (i = 8; i >= 0; i--)
    {
        for (j = 0; j < 30; j++)
        {
            mult_fp(x, x, x1);//squaring for each bit position
            for (k = 0; k < 10; k++)
            {
                x[k] = x1[k];
                x1[k] = 0;
            }
            if ((p_minus_2[i] >> (29 - j)) & 1)//finding the bit from msb to lsb, if the bit is 1 the multiply
            {
                mult_fp(x, a, x1);
                for (k = 0; k < 10; k++)
                {
                    x[k] = x1[k];
                    x1[k] = 0;
                }
            }
        }
    }
    for (k = 0; k < 9; k++)
        a_inv[k] = x[k];
}

void pt_add(long long int *x1, long long int *y1, long long int *x2, long long int *y2, long long int *p_plus_q_x, long long int *p_plus_q_y)
{//lambda = (y_2 - y_1)/(x_2 - x_1), x_3 = lambda^2 - x_1 - x_2, y_3 = lambda * (x_1 - x_3 ) - y1
    long long int p[10] = {1073741823, 1073741823, 1073741823, 63, 0, 0, 4096, 1073725440, 65535, 0};
    long long int subt_y2y1[11] = {0}, subt_y1y2[11] = {0};
    long long int subt_x2x1[11] = {0}, subt_x1x2[11] = {0}, subt_x2x1_inv[11] = {0};
    long long int lambda[10] = {0}, lambda_sq[10] = {0};
    long long int add_x1x2[11] = {0};
    long long int subt_add_x1x2_lambda_sq[10] = {0};
    long long int subt_x1_p_plus_q_x[10] = {0}, subt_p_plus_q_x_x1[10] = {0};
    long long int mult_fp_lambda_subt_x1_p_plus_q_x[10] = {0}, subt_y1_mult_fp_lambda_subt_x1_p_plus_q_x[11] = {0};
    int l;
    if (compare(y2, y1) == 1)
        subt(y2, y1, subt_y2y1, 9);
    else
    {
        subt(y1, y2, subt_y1y2, 9);
        subt(p, subt_y1y2, subt_y2y1, 10);//(y_2 - y_1)
    }
    if (compare(x2, x1) == 1)
        subt(x2, x1, subt_x2x1, 9);//(x_2 - x_1)
    else
    {
        subt(x1, x2, subt_x1x2, 9);
        subt(p, subt_x1x2, subt_x2x1, 10);
    }
    inverse(subt_x2x1, subt_x2x1_inv);//1/(x_2 - x_1)
    mult_fp(subt_y2y1, subt_x2x1_inv, lambda);//(y_2 - y_1)/(x_2 - x_1)
    mult_fp(lambda, lambda, lambda_sq);//lambda^2
    add_fp(x1, x2, add_x1x2);//(x_1 + x_2)
    if (compare(lambda_sq, add_x1x2) == 1)
        subt(lambda_sq, add_x1x2, p_plus_q_x, 9);//lambda^2 - (x_1 + x_2)
    else
    {
        subt(add_x1x2, lambda_sq, subt_add_x1x2_lambda_sq, 9);
        subt(p, subt_add_x1x2_lambda_sq, p_plus_q_x, 10);
    }
    if (compare(x1, p_plus_q_x) == 1)
        subt(x1, p_plus_q_x, subt_x1_p_plus_q_x, 9);//(x_1 - x_3)
    else
    {
        subt(p_plus_q_x, x1, subt_p_plus_q_x_x1, 9);
        subt(p, subt_p_plus_q_x_x1, subt_x1_p_plus_q_x, 10);
    }
    mult_fp(lambda, subt_x1_p_plus_q_x, mult_fp_lambda_subt_x1_p_plus_q_x);//lambda * (x_1 - x_3)
    if (compare(mult_fp_lambda_subt_x1_p_plus_q_x, y1) == 1)
        subt(mult_fp_lambda_subt_x1_p_plus_q_x, y1, p_plus_q_y, 9);//lambda * (x_1 - x_3) - y_1
    else
    {
        subt(y1, mult_fp_lambda_subt_x1_p_plus_q_x, subt_y1_mult_fp_lambda_subt_x1_p_plus_q_x, 9);
        subt(p, subt_y1_mult_fp_lambda_subt_x1_p_plus_q_x, p_plus_q_y, 10);
    }
}

void pt_double(long long int *x1, long long int *y1, long long int *p_plus_q_x, long long int *p_plus_q_y)
{//lambda = (3 * x_1^2 + A) / (2 * y_1)
    long long int p[10] = {1073741823, 1073741823, 1073741823, 63, 0, 0, 4096, 1073725440, 65535, 0};
    long long int a[10] = {1073741820, 1073741823, 1073741823, 63, 0, 0, 4096, 1073725440, 65535, 0};
    long long int add_x1_x1[10] = {0};
    long long int x1_sq[10] = {0};
    long long int three[10] = {3}, mult_fp_three_x1_sq[10] = {0};
    long long int add_mult_fp_three_x1_sq_a[10] = {0};
    long long int two_y1[10] = {0}, two_y1_inv[10] = {0};
    long long int lambda[10] = {0}, lambda_sq[10] = {0};
    long long int add_x1x2[11] = {0};
    long long int subt_add_x1x2_lambda_sq[10] = {0};
    long long int subt_x1_p_plus_q_x[10] = {0}, subt_p_plus_q_x_x1[10] = {0};
    long long int mult_fp_lambda_subt_x1_p_plus_q_x[10] = {0}, subt_y1_mult_fp_lambda_subt_x1_p_plus_q_x[11] = {0};
    int l;
    add_fp(x1, x1, add_x1_x1);
    mult_fp(x1, x1, x1_sq);
    mult_fp(three, x1_sq, mult_fp_three_x1_sq);
    add_fp(mult_fp_three_x1_sq, a, add_mult_fp_three_x1_sq_a);
    add_fp(y1, y1, two_y1);
    inverse(two_y1, two_y1_inv);
    mult_fp(add_mult_fp_three_x1_sq_a, two_y1_inv, lambda);
    mult_fp(lambda, lambda, lambda_sq);//lambda^2
    add_fp(x1, x1, add_x1x2);//(x_1 + x_1)
    if (compare(lambda_sq, add_x1x2) == 1)
        subt(lambda_sq, add_x1x2, p_plus_q_x, 9);//lambda^2 - (x_1 + x_2)
    else
    {
        subt(add_x1x2, lambda_sq, subt_add_x1x2_lambda_sq, 9);
        subt(p, subt_add_x1x2_lambda_sq, p_plus_q_x, 10);
    }
    if (compare(x1, p_plus_q_x) == 1)
        subt(x1, p_plus_q_x, subt_x1_p_plus_q_x, 9);//(x_1 - x_3)
    else
    {
        subt(p_plus_q_x, x1, subt_p_plus_q_x_x1, 9);
        subt(p, subt_p_plus_q_x_x1, subt_x1_p_plus_q_x, 10);
    }
    mult_fp(lambda, subt_x1_p_plus_q_x, mult_fp_lambda_subt_x1_p_plus_q_x);//lambda * (x_1 - x_3)
    if (compare(mult_fp_lambda_subt_x1_p_plus_q_x, y1) == 1)
        subt(mult_fp_lambda_subt_x1_p_plus_q_x, y1, p_plus_q_y, 9);//lambda * (x_1 - x_3) - y_1
    else
    {
        subt(y1, mult_fp_lambda_subt_x1_p_plus_q_x, subt_y1_mult_fp_lambda_subt_x1_p_plus_q_x, 9);
        subt(p, subt_y1_mult_fp_lambda_subt_x1_p_plus_q_x, p_plus_q_y, 10);
    }
}

void scalar_mult(long long int *G_x, long long int *G_y, long long int *A_x, long long int *A_y, long long int *b)
{
    int i, j, k, l, a[2];
    long long int x1[20] = {0}, B_x[10] = {0}, B_y[10] = {0};
    for (i = 0; i < 10; i++)
    {
        A_x[i] = G_x[i];
        A_y[i] = G_y[i];
    }
    first_1(b, a);
    int tk = a[1] + 1;//since first we have to ignore the first one as we are initializing (A_x, A_y) with (G_x, G_y)
    for (i = a[0]; i >= 0; i--)
    {
        for (j = tk; j < 30; j++)
        {
            pt_double(A_x, A_y, B_x, B_y);
            for (k = 0; k < 10; k++)
            {
                A_x[k] = B_x[k];
                A_y[k] = B_y[k];
                B_x[k] = 0;
                B_y[k] = 0;
            }
            if ((b[i] >> (29 - j)) & 1)//if the j-th bit is one then only we have to add
            {
                pt_add(A_x, A_y, G_x, G_y, B_x, B_y);
                for (k = 0; k < 10; k++)
                {
                    A_x[k] = B_x[k];
                    A_y[k] = B_y[k];
                    B_x[k] = 0;
                    B_y[k] = 0;
                }
            }
        }
        tk = 0;//since only for the last limb we have to consider first one
    }
}

int compare(long long int *a, long long int *b) // return 0 if a=b, 1 if a>b, -1 if a<b
{
    int i = 9;
    while (a[i] == b[i] && i >= 0)
        i--;
    if (i == -1)
        return 0;
    if (a[i] > b[i])
        return 1;
    return -1;
}

int first_1(long long int *power, int *a)
{
    int i, j;
    for (i = 8; i >= 0; i--)
    {
        for (j = 0; j < 30; j++)
        {
            if ((power[i] >> (29 - j)) & 1)
            {
                a[0] = i;//value of the limb number of power where first 'one' occur in the binary representation of power
                a[1] = j;//value of the bit number of the limb where first 'one' occur in the binary representation of power
                return 0;
            }
        }
    }
}
