#include "tommath.h"
#include "curve.h"


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
    mp_mulmod(&buf1, &Param->theta, &Param->p, &buf1);   // 3 * theta in buf1
    mp_init_set(&buf2, 4);                               // 4 in buf2
    mp_invmod(&buf2, &Param->p, &buf2);                  // 1/4 in buf2
    mp_mulmod(&buf1, &buf2, &Param->p, &JQ->d);          // 3 * theta / 4 in JQ->d


    /*  Считаем параметр e = -(3*theta^2 + 4 * a) / 16  */
    mp_mulmod(&Param->theta, &Param->theta, &Param->p, &buf1);  // theta^2 in buf1
    mp_clear(&buf2);                                            // clear buf2
    mp_init_set(&buf2, 3);                                      // 3 in buf2
    mp_mulmod(&buf2, &buf1, &Param->p, &buf1);                  // 3 * theta^2 in buf1

    mp_clear(&buf2);                                            // clear buf2
    mp_init_set(&buf2, 4);                                      // 4 in buf2
    mp_mulmod(&buf2, &Param->a, &Param->p, &buf2);              // 4 * a in buf2
    mp_addmod(&buf1, &buf2, &Param->p, &buf1);                  // 3 * theta^2 + 4 * a in buf1
    mp_neg(&buf1, &buf1);                                       // -(3 * theta^2 + 4 * a) in buf1
    mp_add(&buf1, &Param->p, &buf1);                            // -(3 * theta^2 + 4 * a) mod p in buf1
    mp_clear(&buf2);                                            // clear buf2
    mp_init_set(&buf2, 16);                                     // 16 in buf2
    mp_invmod(&buf2, &Param->p, &buf2);                         // 1/16 in buf2
    mp_mulmod(&buf1, &buf2, &Param->p, &JQ->e);                 // -(3 * theta^2 + 4 * a) / 16 in JQ->e


    /*      Считаем координаты порождающего элемента в проективных координатах      */
    /* Координата X = 2 * (x_base - theta) */
    mp_clear(&buf1);                                                // clear buf1
    mp_init_set(&buf1, 2);                                          // 2 in buf1
    mp_submod(&Param->x_base, &Param->theta, &Param->p, &buf2);     // (x_base - theta) in buf2
    mp_mulmod(&buf1, &buf2, &Param->p, &JQ->X);                     // 2 * (x_base - theta) in JQ->X


    /* Координата Y = (2 * x_base + theta) * (x_base - theta)^2 - y_base^2 */
    mp_clear(&buf1);                                               // clear buf1
    mp_init_set(&buf1, 2);                                         // 2 in buf1
    mp_mulmod(&buf1, &Param->x_base, &Param->p, &buf1);            // 2 * x_base in buf1
    mp_addmod(&buf1, &Param->theta, &Param->p, &buf1);             // 2 * x_base + theta in buf1

    mp_submod(&Param->x_base, &Param->theta, &Param->p, &buf2);    // x_base - theta in buf2
    mp_mulmod(&buf2, &buf2, &Param->p, &buf2);                     // (x_base - theta)^2 in buf2
    mp_mulmod(&buf1, &buf2, &Param->p, &buf1);                     // (2 * x_base + theta) * (x_base - theta)^2 in buf1

    mp_mulmod(&Param->y_base, &Param->y_base, &Param->p, &buf2);   // y_base^2 in buf2
    mp_submod(&buf1, &buf2, &Param->p, &JQ->Y);                    // (2 * x_base + theta) * (x_base - theta)^2 - y_base^2 in JQ->Y

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
    int out_size = 250;
    char out[out_size];

    mp_to_radix(&P->X, out, out_size, NULL, 10);
    printf("X = %s\n", out);
    mp_to_radix(&P->Y, out, out_size, NULL, 10);
    printf("Y = %s\n", out);
    mp_to_radix(&P->Z, out, out_size, NULL, 10);
    printf("Z = %s\n", out);
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
    mp_mulmod(&P->X, &buf, &JQ->p, &x);     // X / Z in x

    /* y = Y / Z^2 */
    mp_mulmod(&buf, &buf, &JQ->p, &buf);    // 1 / Z^2 in buf
    mp_mulmod(&P->Y, &buf, &JQ->p, &y);     // Y / Z^2 in y

    printf("Точка в афинных координатах:\n");
    int out_size = 250;
    char out[out_size];

    mp_to_radix(&x, out, out_size, NULL, 10);
    printf("x = %s\n", out);
    mp_to_radix(&y, out, out_size, NULL, 10);
    printf("y = %s\n", out);

    mp_clear_multi(&x, &y, &buf, NULL);
}


/*
 * Основные функции
 */
// Проверить находится ли данная точка на данной кривой
int IsPointOnCurve(struct Point * P, struct JacobiQuadric * JQ)
{
    mp_int buf1, buf2, buf3;
    mp_init_multi(&buf2, &buf3, NULL);
    int IsOnCurve;

    /* Для проверки нахождения точки на кривой подставляем её коориднаты (X : Y: Z)
     * в уравение кривой Y ^ 2 == e * X^4 - 2 * d * X^2 * Y^2 + Z^4              */
    mp_init_set(&buf1, 2);                      // 2 in buf1
    mp_exptmod(&P->X, &buf1, &JQ->p, &buf2);    // X^2 in buf2
    mp_exptmod(&P->Z, &buf1, &JQ->p, &buf3);    // Z^2 in buf3
    mp_mulmod(&buf2, &buf3, &JQ->p, &buf3);     // X^2 * Z^2 in buf3
    mp_mulmod(&buf1, &JQ->d, &JQ->p, &buf1);    // 2 * d in buf1
    mp_mulmod(&buf1, &buf3, &JQ->p, &buf3);     // X^2 * Z^2 * 2 * d in buf3

    mp_clear(&buf1);                            // clear buf1
    mp_init_set(&buf1, 4);                      // 4 in buf1
    mp_exptmod(&P->X, &buf1, &JQ->p, &buf2);    // X^4 in buf2
    mp_mulmod(&buf2, &JQ->e, &JQ->p, &buf2);    // e * X^4 in buf2

    mp_submod(&buf2, &buf3, &JQ->p, &buf3);     // e * X^4 - X^2 * Y^2 * 2 * d in buf3

    mp_exptmod(&P->Z, &buf1, &JQ->p, &buf2);    // Z^4 in buf2
    mp_addmod(&buf3, &buf2, &JQ->p, &buf3);     // e * X^4 - 2 * d * X^2 * Y^2 + Z^4 in buf3

    mp_clear(&buf1);                            // clear buf1
    mp_init_set(&buf1, 2);                      // 2 in buf1
    mp_exptmod(&P->Y, &buf1, &JQ->p, &buf2);    // Y^2 in buf2

    if (mp_cmp(&buf2, &buf3) == MP_EQ)
    {
        IsOnCurve = 0;
    }
    else
    {
        IsOnCurve = -1;
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

    mp_mulmod(&T1, &T3, &JQ->p, &T7);       // T7 = T1 * T3
    mp_addmod(&T2, &T7, &JQ->p, &T7);       // T7 = T2 + T7
    mp_mulmod(&T4, &T6, &JQ->p, &T8);       // T8 = T4 * T6
    mp_addmod(&T5, &T8, &JQ->p, &T8);       // T8 = T5 + T8
    mp_mulmod(&T2, &T5, &JQ->p, &T2);       // T2 = T2 * T5
    mp_mulmod(&T7, &T8, &JQ->p, &T7);       // T7 = T7 * T8
    mp_submod(&T7, &T2, &JQ->p, &T7);       // T7 = T7 - T2
    mp_mulmod(&T1, &T4, &JQ->p, &T5);       // T5 = T1 * T4
    mp_addmod(&T1, &T3, &JQ->p, &T1);       // T1 = T1 + T3
    mp_mulmod(&T3, &T6, &JQ->p, &T8);       // T8 = T3 * T6
    mp_addmod(&T4, &T6, &JQ->p, &T4);       // T4 = T4 + T6
    mp_mulmod(&T5, &T8, &JQ->p, &T6);       // T6 = T5 * T8
    mp_submod(&T7, &T6, &JQ->p, &T7);       // T7 = T7 - T6
    mp_mulmod(&T1, &T4, &JQ->p, &T1);       // T1 = T1 * T4
    mp_submod(&T1, &T5, &JQ->p, &T1);       // T1 = T1 - T5
    mp_submod(&T1, &T8, &JQ->p, &T1);       // T1 = T1 - T8
    mp_mulmod(&T1, &T1, &JQ->p, &T3);       // T3 = T1 * T1
    mp_addmod(&T6, &T6, &JQ->p, &T6);       // T6 = T6 + T6
    mp_submod(&T3, &T6, &JQ->p, &T3);       // T3 = T3 - T6
    mp_mulmod(&JQ->e, &T6, &JQ->p, &T4);    // T4 = e * T6
    mp_mulmod(&T3, &T4, &JQ->p, &T3);       // T3 = T3 * T4
    mp_mulmod(&JQ->d, &T6, &JQ->p, &T4);    // T4 = d * T6
    mp_submod(&T2, &T4, &JQ->p, &T2);       // T2 = T2 - T4
    mp_mulmod(&T8, &T8, &JQ->p, &T4);       // T4 = T8 * T8
    mp_mulmod(&T5, &T5, &JQ->p, &T8);       // T8 = T5 * T5
    mp_mulmod(&JQ->e, &T8, &JQ->p, &T8);    // T8 = e * T8
    mp_addmod(&T4, &T8, &JQ->p, &T5);       // T5 = T4 + T8
    mp_mulmod(&T2, &T5, &JQ->p, &T2);       // T2 = T2 * T5
    mp_addmod(&T2, &T3, &JQ->p, &T2);       // T2 = T2 + T3
    mp_submod(&T4, &T8, &JQ->p, &T5);       // T5 = T4 - T8

    mp_clear_multi(&P3->X, &P3->Y, &P3->Z, NULL);
    mp_init_copy(&P3->X, &T7);              // X3 = T7
    mp_init_copy(&P3->Y, &T2);              // Y3 = T2
    mp_init_copy(&P3->Z, &T5);              //Z3 = T5

    mp_clear_multi(&T1, &T2, &T3, &T4, &T5, &T6, &T7, &T8, NULL);
}


// Проверка равенства двух точек на кривой
int ArePointsEqual(struct Point * P1, struct Point * P2, struct JacobiQuadric * JQ)
{
    mp_int x1, y1, x2, y2, buf;
    mp_init_multi(&x1, &y1, &x2, &y2, &buf, NULL);

    /*  Перевод первой точки в афинные координаты   */
    /* x = X / Z */
    mp_invmod(&P1->Z, &JQ->p, &buf);          // 1 / Z in buf
    mp_mulmod(&P1->X, &buf, &JQ->p, &x1);     // X / Z in x1

    /* y = Y / Z^2 */
    mp_mulmod(&buf, &buf, &JQ->p, &buf);      // 1 / Z^2 in buf
    mp_mulmod(&P1->Y, &buf, &JQ->p, &y1);     // Y / Z^2 in y1

    /*  Перевод второй точки в афинные координаты   */
    /* x = X / Z */
    mp_invmod(&P2->Z, &JQ->p, &buf);          // 1 / Z in buf
    mp_mulmod(&P2->X, &buf, &JQ->p, &x2);     // X / Z in x2

    /* y = Y / Z^2 */
    mp_mulmod(&buf, &buf, &JQ->p, &buf);     // 1 / Z^2 in buf
    mp_mulmod(&P2->Y, &buf, &JQ->p, &y2);    // Y / Z^2 in y2

    /*  Сравнение точек в афинных координатах   */
    if ((mp_cmp(&x1, &x2) == MP_EQ) && (mp_cmp(&y1, &y2) == MP_EQ))
    {
        mp_clear_multi(&x1, &y1, &x2, &y2, &buf, NULL);
        return 0;
    }
    else
    {
        mp_clear_multi(&x1, &y1, &x2, &y2, &buf, NULL);
        return -1;
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

    int n = mp_count_bits(k);

    // Лесенка Монтгомери
    int i = 0;
    for (i = 0; i < n; i++)
    {
        if (mp_get_bit(k, n - i - 1) == 1)
        {
            Addition(&R, &Q, &Q, JQ);   // Q = R + Q
            Addition(&R, &R, &R, JQ);   // R = R + R
        }
        else
        {
            Addition(&R, &Q, &R, JQ);   // R = R + Q
            Addition(&Q, &Q, &Q, JQ);   // Q = Q + Q
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

