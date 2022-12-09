#include "EC_p256.c"

void main()
{
    clock_t start = clock();
    long long int gx[10] = {412664470, 310699287, 515062287, 14639179, 608236151, 865834382, 69500811, 880588875, 27415};//x-co-ordinate of the generator
    long long int gy[10] = {935285237, 785973664, 857074924, 864867802, 262018603, 531442160, 670677230, 280543110, 20451};//y-co-ordinate of the generator
    long long int n[10] = {0x3c632551, 0xee72b0b, 0x3179e84f, 0x39beab69, 0x3fffffbc, 0x3fffffff, 0xfff, 0x3fffc000, 0xffff};//this is the 2^30 base representation of the order of the ellitic curve
    long long int a_gx[10] = {0}, a_gy[10] = {0};
    long long int b_gx[10] = {0}, b_gy[10] = {0};
    long long int ab_gx[10] = {0}, ab_gy[10] = {0};
    long long int ba_gx[10] = {0}, ba_gy[10] = {0};
    long long int a[10] = {0}, b[10] = {0};
    int i, l;
    srand(time(NULL));
    for(i = 0; i < 9; i++)//assigning random values to the scalar, at the same time making sure the value of the scalar must be less than the order of the elliptic curve
    {
        a[i] = rand() & n[i];
        b[i] = rand() & n[i];
    }
    printf("a: ");
    for (l = 0; l < 9; l++)
        printf("%10lld\t", a[l]);
    printf("\n");
    printf("b: ");
    for (l = 0; l < 9; l++)
        printf("%10lld\t", b[l]);
    printf("\n\n");
    scalar_mult(gx, gy, a_gx, a_gy, a);//aG
    printf("a_gx: ");
    for (l = 0; l < 9; l++)
        printf("%10lld\t", a_gx[l]);
    printf("\n");
    printf("a_gy: ");
    for (l = 0; l < 9; l++)
        printf("%10lld\t", a_gy[l]);
    printf("\n\n");
    scalar_mult(a_gx, a_gy, ba_gx, ba_gy, b);
    printf("ba_gx: ");//baG
    for (l = 0; l < 9; l++)
        printf("%10lld\t", ba_gx[l]);
    printf("\n");
    printf("ba_gy: ");
    for (l = 0; l < 9; l++)
        printf("%10lld\t", ba_gy[l]);
    printf("\n\n");
    scalar_mult(gx, gy, b_gx, b_gy, b);//bG
    printf("b_gx: ");
    for (l = 0; l < 9; l++)
        printf("%10lld\t", b_gx[l]);
    printf("\n");
    printf("b_gy: ");
    for (l = 0; l < 9; l++)
        printf("%10lld\t", b_gy[l]);
    printf("\n\n");
    scalar_mult(b_gx, b_gy, ab_gx, ab_gy, a);//abG
    printf("ab_gx: ");
    for (l = 0; l < 9; l++)
        printf("%10lld\t", ab_gx[l]);
    printf("\n");
    printf("ab_gy: ");
    for (l = 0; l < 9; l++)
        printf("%10lld\t", ab_gy[l]);
    printf("\n\n");
    clock_t end = clock();
    double elapsed = ((float)end - (float)start) / CLOCKS_PER_SEC;
    printf("Time measured: %.3f seconds.\n", elapsed);
}