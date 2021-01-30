#include "tommath.h"
#include "curve.h"

/*
 * Дополнительные функции, обеспечивающие арифметику по модулю
 */
// Находит сумму а + b по модулю mod и помещает результат в с
void mp_add_mod(mp_int * a, mp_int * b, mp_int * c, mp_int * mod)
{
    mp_add(a, b, c);
    mp_div(c, mod, NULL, c);
}

// Находит разность а - b по модулю mod и помещает результат в с
void mp_sub_mod(mp_int * a, mp_int * b, mp_int * c, mp_int * mod)
{
    mp_sub(a, b, c);
    mp_int ZERO;
    mp_init_set(&ZERO, 0);

    if (mp_cmp(c, &ZERO) == MP_GT)
    {
        // c > 0
        mp_div(c, mod, NULL, c);
    }
    else
    {
        // c <= 0
        mp_add(c, mod, c);
        mp_div(c, mod, NULL, c);
    }

    mp_clear(&ZERO);
}

// Находит произведение а * b по модулю mod и помещает результат в с
void mp_mul_mod(mp_int * a, mp_int * b, mp_int * c, mp_int * mod)
{
    mp_mul(a, b, c);
    mp_div(c, mod, NULL, c);
}

// Считает -а по модулю mod и помещает результат в b
void mp_neg_mod(mp_int * a, mp_int * b, mp_int * mod)
{
    mp_neg(a, b);
    mp_add(b, mod, b);
}


/*
 * Функции инициализации основных структур
 */
// Инициализация точки
void InitPoint(struct Point * P, mp_int * x, mp_int * y, mp_int * z)
{
    mp_init_copy(&P->X, x);
    mp_init_copy(&P->Y, y);
    mp_init_copy(&P->Z, z);
}

// Инициализация структуры с параметрами
void InitParameters(struct Parameters * Param, mp_int * p, mp_int * q, mp_int * a, mp_int * b, mp_int * x_base, mp_int * y_base, mp_int * theta)
{
    mp_init_copy(&Param->p, p);
    mp_init_copy(&Param->q, q);

    mp_init_copy(&Param->a, a);
    mp_init_copy(&Param->b, b);

    mp_init_copy(&Param->x_base, x_base);
    mp_init_copy(&Param->y_base, y_base);

    mp_init_copy(&Param->theta, theta);

}

// Инициализация структуры с параметрами кривой
void InitJacobiQuadric(struct JacobiQuadric * JQ, struct Parameters * Param)
{
    mp_init_copy(&JQ->p, &Param->p);
    mp_init_multi(&JQ->d, &JQ->e, &JQ->X, &JQ->Y, NULL);


    /*      Считаем параметры e и d квадрики Якоби       */
    /*  Считаем параметр d = 3 * theta / 4  */
    mp_int buf1, buf2;

    mp_init_set(&buf1, 3);                               // 3 in buf1
    mp_mul_mod(&buf1, &Param->theta, &buf1, &Param->p);  // 3 * theta in buf1
    mp_init_set(&buf2, 4);                               // 4 in buf2
    mp_invmod(&buf2, &Param->p, &buf2);                  // 1/4 in buf2
    mp_mul_mod(&buf1, &buf2, &JQ->d, &Param->p);         // 3 * theta / 4 in JQ->d


    /*  Считаем параметр e = -(3*theta^2 + 4 * a) / 16  */
    mp_mul_mod(&Param->theta, &Param->theta, &buf1, &Param->p); // theta^2 in buf1
    mp_clear(&buf2);                                            // clear buf2
    mp_init_set(&buf2, 3);                                      // 3 in buf2
    mp_mul_mod(&buf2, &buf1, &buf1, &Param->p);                 // 3 * theta^2 in buf1

    mp_clear(&buf2);                                            // clear buf2
    mp_init_set(&buf2, 4);                                      // 4 in buf2
    mp_mul_mod(&buf2, &Param->a, &buf2, &Param->p);             // 4 * a in buf2
    mp_add_mod(&buf1, &buf2, &buf1, &Param->p);                 // 3 * theta^2 + 4 * a in buf1

    mp_neg_mod(&buf1, &buf1, &Param->p);                        // -(3 * theta^2 + 4 * a) in buf1
    mp_clear(&buf2);                                            // clear buf2
    mp_init_set(&buf2, 16);                                     // 16 in buf2
    mp_invmod(&buf2, &Param->p, &buf2);                         // 1/16 in buf2
    mp_mul_mod(&buf1, &buf2, &JQ->e, &Param->p);                // -(3 * theta^2 + 4 * a) / 16 in JQ->e


    /*      Считаем координаты порождающего элемента в проективных координатах      */
    /* Координата X = 2 * (x_base - theta) */
    mp_clear(&buf1);                                                // clear buf1
    mp_init_set(&buf1, 2);                                          // 2 in buf1
    mp_sub_mod(&Param->x_base, &Param->theta, &buf2, &Param->p);    // (x_base - theta) in buf2
    mp_mul_mod(&buf1, &buf2, &JQ->X, &Param->p);                    // 2 * (x_base - theta) in JQ->X


    /* Координата Y = (2 * x_base + theta) * (x_base - theta)^2 - y_base^2 */
    mp_clear(&buf1);                                               // clear buf1
    mp_init_set(&buf1, 2);                                         // 2 in buf1
    mp_mul_mod(&buf1, &Param->x_base, &buf1, &Param->p);           // 2 * x_base in buf1
    mp_add_mod(&buf1, &Param->theta, &buf1, &Param->p);            // 2 * x_base + theta in buf1

    mp_sub_mod(&Param->x_base, &Param->theta, &buf2, &Param->p);   // x_base - theta in buf2
    mp_mul_mod(&buf2, &buf2, &buf2, &Param->p);                    // (x_base - theta)^2 in buf2
    mp_mul_mod(&buf1, &buf2, &buf1, &Param->p);                    // (2 * x_base + theta) * (x_base - theta)^2 in buf1

    mp_mul_mod(&Param->y_base, &Param->y_base, &buf2, &Param->p);  // y_base^2 in buf2
    mp_sub_mod(&buf1, &buf2, &JQ->Y, &Param->p);                   // (2 * x_base + theta) * (x_base - theta)^2 - y_base^2 in JQ->Y

    mp_clear_multi(&buf1, &buf2, NULL);                            // очистка переменных - буферов


    /* Координата Z совпадает с координатой y в афинной форме: Z = y_base */
    mp_init_copy(&JQ->Z, &Param->y_base);
}


/*
 * Функции освобождения памяти, использующейся для хранения основных структур
 */
// Очистка памяти, использовавшейся для хранения точки
void ClearPoint(struct Point * P)
{
    mp_clear_multi(&P->X, &P->Y, &P->Z, NULL);
}

// Очистка памяти, использовавшейся для хранения параметров
void ClearParameters(struct Parameters * Param)
{
    mp_clear_multi(&Param->p, &Param->q, &Param->a, &Param->b, &Param->x_base, &Param->y_base, &Param->theta, NULL);
}

// Очистка памяти, использовавшейся для хранения параметров кривой Якоби
void ClearJacobiQuadric(struct JacobiQuadric * JQ)
{
    mp_clear_multi(&JQ->p, &JQ->e, &JQ->d, &JQ->X, &JQ->Y, &JQ->Z, NULL);
}


/*
 * Функции для вывода координат точки на экран
 */
// Вывести значения проективных координат точки на экран
void PrintPoint(struct Point * P)
{
    printf("Точка в проективных координатах:\n");
    size_t size_x, size_y, size_z;

    mp_radix_size(&P->X, 10, &size_x);
    mp_radix_size(&P->Y, 10, &size_y);
    mp_radix_size(&P->Z, 10, &size_z);

    char x_str[size_x], y_str[size_y], z_str[size_z];

    mp_to_radix(&P->X, x_str, size_x, NULL, 10);
    printf("X = %s\n", x_str);
    mp_to_radix(&P->Y, y_str, size_y, NULL, 10);
    printf("Y = %s\n", y_str);
    mp_to_radix(&P->Z, z_str, size_z, NULL, 10);
    printf("Z = %s\n", z_str);
}

// Вывести значения афинных координат точки на экран
void PrintPointAffine(struct Point * P, struct JacobiQuadric * JQ)
{
    /* Афинные координаты точки (x : y) выражаются через проективные (X : Y : Z)
     * следующим образом: x = X / Z; y = Y / Z^2        */
    mp_int x, y, buf;
    mp_init_multi(&x, &y, &buf, NULL);

    /* x = X / Z */
    mp_invmod(&P->Z, &JQ->p, &buf);         // 1 / Z in buf
    mp_mul_mod(&P->X, &buf, &x, &JQ->p);    // X / Z in x

    /* y = Y / Z^2 */
    mp_mul_mod(&buf, &buf, &buf, &JQ->p);   // 1 / Z^2 in buf
    mp_mul_mod(&P->Y, &buf, &y, &JQ->p);    // Y / Z^2 in y

    printf("Точка в афинных координатах:\n");
    size_t size_x, size_y;

    mp_radix_size(&x, 10, &size_x);
    mp_radix_size(&y, 10, &size_y);;

    char x_str[size_x], y_str[size_y];

    mp_to_radix(&x, x_str, size_x, NULL, 10);
    printf("x = %s\n", x_str);
    mp_to_radix(&y, y_str, size_y, NULL, 10);
    printf("y = %s\n", y_str);

    mp_clear_multi(&x, &y, &buf, NULL);
}


/*
 * Основные функции
 */
// Проверить находится ли данная точка на данной кривой
bool IsPointOnCurve(struct Point * P, struct JacobiQuadric * JQ)
{
    mp_int buf1, buf2, buf3;
    mp_init_multi(&buf2, &buf3, NULL);
    bool IsOnCurve;

    /* Для проверки нахождения точки на кривой подставляем её коориднаты (X : Y: Z)
     * в уравение кривой Y ^ 2 == e * X^4 - 2 * d * X^2 * Y^2 + Z^4              */
    mp_init_set(&buf1, 2);                      // 2 in buf1
    mp_exptmod(&P->X, &buf1, &JQ->p, &buf2);    // X^2 in buf2
    mp_exptmod(&P->Z, &buf1, &JQ->p, &buf3);    // Z^2 in buf3
    mp_mul_mod(&buf2, &buf3, &buf3, &JQ->p);    // X^2 * Z^2 in buf3
    mp_mul_mod(&buf1, &JQ->d, &buf1, &JQ->p);   // 2 * d in buf1
    mp_mul_mod(&buf1, &buf3, &buf3, &JQ->p);    // X^2 * Z^2 * 2 * d in buf3

    mp_clear(&buf1);                            // clear buf1
    mp_init_set(&buf1, 4);                      // 4 in buf1
    mp_exptmod(&P->X, &buf1, &JQ->p, &buf2);    // X^4 in buf2
    mp_mul_mod(&buf2, &JQ->e, &buf2, &JQ->p);   // e * X^4 in buf2

    mp_sub_mod(&buf2, &buf3, &buf3, &JQ->p);    // e * X^4 - X^2 * Y^2 * 2 * d in buf3

    mp_exptmod(&P->Z, &buf1, &JQ->p, &buf2);    // Z^4 in buf2
    mp_add_mod(&buf3, &buf2, &buf3, &JQ->p);    // e * X^4 - 2 * d * X^2 * Y^2 + Z^4 in buf3

    mp_clear(&buf1);                            // clear buf1
    mp_init_set(&buf1, 2);                      // 2 in buf1
    mp_exptmod(&P->Y, &buf1, &JQ->p, &buf2);    // Y^2 in buf2

    if (mp_cmp(&buf2, &buf3) == MP_EQ)
    {
        IsOnCurve = true;
        printf("Точка принадлежит кривой\n");
    }
    else
    {
        IsOnCurve = false;
        printf("Точка не принадлежит кривой\n");
    }

    mp_clear_multi(&buf1, &buf2, &buf3, NULL);

    return IsOnCurve;
}


// Сложение двух точек P1 + P2, результат записывается в третью точку P3
void Addition(struct Point * P1, struct Point * P2, struct Point * P3, struct JacobiQuadric * JQ)
{
    mp_int T1, T2, T3, T4, T5, T6, T7, T8;

    mp_init_multi(&T7, &T8, NULL);
    mp_init_copy(&T1, &P1->X);              // T1 = X1
    mp_init_copy(&T2, &P1->Y);              // T2 = Y1
    mp_init_copy(&T3, &P1->Z);              // T3 = Z1
    mp_init_copy(&T4, &P2->X);              // T4 = X2
    mp_init_copy(&T5, &P2->Y);              // T5 = Y2
    mp_init_copy(&T6, &P2->Z);              // T6 = Z6

    mp_mul_mod(&T1, &T3, &T7, &JQ->p);      // T7 = T1 * T3
    mp_add_mod(&T2, &T7, &T7, &JQ->p);      // T7 = T2 + T7
    mp_mul_mod(&T4, &T6, &T8, &JQ->p);      // T8 = T4 * T6
    mp_add_mod(&T5, &T8, &T8, &JQ->p);      // T8 = T5 + T8
    mp_mul_mod(&T2, &T5, &T2, &JQ->p);      // T2 = T2 * T5
    mp_mul_mod(&T7, &T8, &T7, &JQ->p);      // T7 = T7 * T8
    mp_sub_mod(&T7, &T2, &T7, &JQ->p);      // T7 = T7 - T2
    mp_mul_mod(&T1, &T4, &T5, &JQ->p);      // T5 = T1 * T4
    mp_add_mod(&T1, &T3, &T1, &JQ->p);      // T1 = T1 + T3
    mp_mul_mod(&T3, &T6, &T8, &JQ->p);      // T8 = T3 * T6
    mp_add_mod(&T4, &T6, &T4, &JQ->p);      // T4 = T4 + T6
    mp_mul_mod(&T5, &T8, &T6, &JQ->p);      // T6 = T5 * T8
    mp_sub_mod(&T7, &T6, &T7, &JQ->p);      // T7 = T7 - T6
    mp_mul_mod(&T1, &T4, &T1, &JQ->p);      // T1 = T1 * T4
    mp_sub_mod(&T1, &T5, &T1, &JQ->p);      // T1 = T1 - T5
    mp_sub_mod(&T1, &T8, &T1, &JQ->p);      // T1 = T1 - T8
    mp_mul_mod(&T1, &T1, &T3, &JQ->p);      // T3 = T1 * T1
    mp_add_mod(&T6, &T6, &T6, &JQ->p);      // T6 = T6 + T6
    mp_sub_mod(&T3, &T6, &T3, &JQ->p);      // T3 = T3 - T6
    mp_mul_mod(&JQ->e, &T6, &T4, &JQ->p);   // T4 = e * T6
    mp_mul_mod(&T3, &T4, &T3, &JQ->p);      // T3 = T3 * T4
    mp_mul_mod(&JQ->d, &T6, &T4, &JQ->p);   // T4 = d * T6
    mp_sub_mod(&T2, &T4, &T2, &JQ->p);      // T2 = T2 - T4
    mp_mul_mod(&T8, &T8, &T4, &JQ->p);      // T4 = T8 * T8
    mp_mul_mod(&T5, &T5, &T8, &JQ->p);      // T8 = T5 * T5
    mp_mul_mod(&JQ->e, &T8, &T8, &JQ->p);   // T8 = e * T8
    mp_add_mod(&T4, &T8, &T5, &JQ->p);      // T5 = T4 + T8
    mp_mul_mod(&T2, &T5, &T2, &JQ->p);      // T2 = T2 * T5
    mp_add_mod(&T2, &T3, &T2, &JQ->p);      // T2 = T2 + T3
    mp_sub_mod(&T4, &T8, &T5, &JQ->p);      // T5 = T4 - T8

    mp_clear_multi(&P3->X, &P3->Y, &P3->Z, NULL);
    mp_init_copy(&P3->X, &T7);              // X3 = T7
    mp_init_copy(&P3->Y, &T2);              // Y3 = T2
    mp_init_copy(&P3->Z, &T5);              //Z3 = T5

    mp_clear_multi(&T1, &T2, &T3, &T4, &T5, &T6, &T7, &T8, NULL);
}


// Проверка равенства двух точек на кривой
bool ArePointsEqual(struct Point * P1, struct Point * P2, struct JacobiQuadric * JQ)
{
    mp_int x1, y1, x2, y2, buf;
    mp_init_multi(&x1, &y1, &x2, &y2, &buf, NULL);

    /*  Перевод первой точки в афинные координаты   */
    /* x = X / Z */
    mp_invmod(&P1->Z, &JQ->p, &buf);          // 1 / Z in buf
    mp_mul_mod(&P1->X, &buf, &x1, &JQ->p);    // X / Z in x1

    /* y = Y / Z^2 */
    mp_mul_mod(&buf, &buf, &buf, &JQ->p);     // 1 / Z^2 in buf
    mp_mul_mod(&P1->Y, &buf, &y1, &JQ->p);    // Y / Z^2 in y1

    /*  Перевод второй точки в афинные координаты   */
    /* x = X / Z */
    mp_invmod(&P2->Z, &JQ->p, &buf);          // 1 / Z in buf
    mp_mul_mod(&P2->X, &buf, &x2, &JQ->p);    // X / Z in x2

    /* y = Y / Z^2 */
    mp_mul_mod(&buf, &buf, &buf, &JQ->p);     // 1 / Z^2 in buf
    mp_mul_mod(&P2->Y, &buf, &y2, &JQ->p);    // Y / Z^2 in y1

    if ((mp_cmp(&x1, &x2) == MP_EQ) && (mp_cmp(&y1, &y2) == MP_EQ))
    {
        mp_clear_multi(&x1, &y1, &x2, &y2, &buf, NULL);
        return true;
    }
    else
    {
        mp_clear_multi(&x1, &y1, &x2, &y2, &buf, NULL);
        return false;
    }
}


// Реализация алгоритма "лесенка Монтгомери"
void MontgomeryLadder(struct Point * P, mp_int * k, struct Point * RES, struct JacobiQuadric * JQ)
{
    // Q = (0 : 1 : 1)
    struct Point Q;
    mp_int qx, qy, qz;
    mp_init_set(&qx, 0);
    mp_init_set(&qy, 1);
    mp_init_set(&qz, 1);
    InitPoint(&Q, &qx, &qy, &qz);

    // R = P
    struct Point R;
    mp_int rx, ry, rz;
    mp_init_copy(&rx, &P->X);
    mp_init_copy(&ry, &P->Y);
    mp_init_copy(&rz, &P->Z);
    InitPoint(&R, &rx, &ry, &rz);

    size_t n;
    mp_radix_size(k, 2, &n);
    char bin_k[n];
    mp_to_radix(k, bin_k, n, NULL, 2);

    // Лесенка Монтгомери
    int i = 0;
    for (i = 0; i < n - 1; i++)
    {
        if (bin_k[i] == '1')
        {
            Addition(&R, &Q, &Q, JQ);   // Q = R + Q
            Addition(&R, &R, &R, JQ);   // R = R + R
        }
        else
        {
            Addition(&R, &Q, &R, JQ); // R = R + Q
            Addition(&Q, &Q, &Q, JQ); // Q = Q + Q
        }
    }

    mp_clear_multi(&RES->X, &RES->Y, &RES->Z, NULL);
    mp_init_copy(&RES->X, &Q.X);
    mp_init_copy(&RES->Y, &Q.Y);
    mp_init_copy(&RES->Z, &Q.Z);

    mp_clear_multi(&rx, &ry, &rz, NULL);
    mp_clear_multi(&qx, &qy, &qz, NULL);
    ClearPoint(&R);
    ClearPoint(&Q);
}

