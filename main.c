#include <stdio.h>
#include "tommath.h"
#include "curve.h"

void Test()
{
    // Инициализация параметров
    mp_int p, q, a, b, x_base, y_base, theta;
    mp_init_multi(&p, &a, &b, &q, &x_base, &y_base, &theta, NULL);
    mp_read_radix(&p, p_str, 10);
    mp_read_radix(&a, a_str, 10);
    mp_read_radix(&b, b_str, 10);
    mp_read_radix(&q, q_str, 10);
    mp_read_radix(&x_base, x_base_str, 10);
    mp_read_radix(&y_base, y_base_str, 10);
    mp_read_radix(&theta, theta_str, 10);

    struct Parameters ParamSet;
    InitParameters(&ParamSet, &p, &q, &a, &b, &x_base, &y_base, &theta);

    struct JacobiQuadric Quadric;
    InitJacobiQuadric(&Quadric, &ParamSet);

    // Вывод параметров
    printf("Посчитанные параметры квадрики Якоби:\n");
    char answer[250];
    mp_to_radix(&Quadric.d, answer, 250, NULL, 10);
    printf("d = %s\n", answer);

    mp_to_radix(&Quadric.e, answer, 250, NULL, 10);
    printf("e = %s\n", answer);

    mp_to_radix(&Quadric.X, answer, 250, NULL, 10);
    printf("X_base = %s\n", answer);

    mp_to_radix(&Quadric.Y, answer, 250, NULL, 10);
    printf("Y_base = %s\n", answer);

    mp_to_radix(&Quadric.Z, answer, 250, NULL, 10);
    printf("Z_base = %s\n", answer);

    // Вводим порождающий элемент (P_base) и нейтральный элемент Е
    struct Point P_base;
    struct Point E;
    mp_int ex, ey, ez;
    mp_init_set(&ex, 0);
    mp_init_set(&ey, 1);
    mp_init_set(&ez, 1);

    InitPoint(&P_base, &Quadric.X, &Quadric.Y, &Quadric.Z);
    InitPoint(&E, &ex, &ey, &ez);


    // ТЕСТИРОВАНИЕ //
    bool ALL_TESTS_ARE_CORRECT = true;

    // ТЕСТ 1: ПРОВЕРКА ПРИНАДЛЕЖНОСТИ НЕЙТРАЛЬНОГО ЭЛЕМЕНТА
    printf("\nТЕСТ 1: ПРОВЕРКА ПРИНАДЛЕЖНОСТИ НЕЙТРАЛЬНОГО ЭЛЕМЕНТА\n");
    PrintPoint(&E);
    PrintPointAffine(&E, &Quadric);
    if (IsPointOnCurve(&E, &Quadric))
    {
        printf("Нейтральный элемент Е принадлежит кривой\n");
    }
    else
    {
        printf("Нейтральный элемент Е не принадлежит кривой\n");
        ALL_TESTS_ARE_CORRECT = false;
    }


    // ТЕСТ 2: ПОРОЖДАЮЩИЙ ЭЛЕМЕНТ В АФИННЫХ КООРДИНАТАХ
    printf("\nТЕСТ 2: ПОРОЖДАЮЩИЙ ЭЛЕМЕНТ В АФИННЫХ КООРДИНАТАХ\n");
    PrintPointAffine(&P_base, &Quadric);
    if (IsPointOnCurve(&P_base, &Quadric))
    {
        printf("Порождающий элемент P_base принадлежит кривой\n");
    }
    else
    {
        printf("Порождающий элемент P_base не принадлежит кривой\n");
        ALL_TESTS_ARE_CORRECT = false;
    }


    // ТECT 3: E + P_base = P_base?
    printf("\nТECT 3: E + P_base = P_base?\n");
    struct Point SumPoint;
    mp_int spx, spy, spz;
    mp_init_multi(&spx, &spy, &spz, NULL);
    InitPoint(&SumPoint, &spx, &spy, &spz);
    Addition(&E, &P_base, &SumPoint, &Quadric);

    PrintPoint(&SumPoint);
    PrintPointAffine(&SumPoint, &Quadric);
    if (IsPointOnCurve(&SumPoint, &Quadric))
    {
        printf("Точка Е + P_base принадлежит кривой\n");
    }
    else
    {
        printf("Точка Е + P_base не принадлежит кривой\n");
        ALL_TESTS_ARE_CORRECT = false;
    }

    if (ArePointsEqual(&SumPoint, &P_base, &Quadric))
    {
        printf("Точки E+P_base и P_base равны\n");
    }
    else
    {
        printf("Точки E+P_base и P_base не равны\n");
        ALL_TESTS_ARE_CORRECT = false;
    }


    // ТECT 4: Принадлежит ли точка P2 = (5 : 1 : 4) кривой
    printf("\nТECT 4: Принадлежит ли точка P2 = (5 : 1 : 4) кривой\n");
    struct Point P2;
    mp_int p2x, p2y, p2z;
    mp_init_set(&p2x, 5);
    mp_init_set(&p2y, 1);
    mp_init_set(&p2z, 4);
    InitPoint(&P2, &p2x, &p2y, &p2z);

    PrintPointAffine(&P2, &Quadric);
    if (IsPointOnCurve(&P2, &Quadric))
    {
        printf("Точка (5 : 1 : 4) принадлежит кривой\n");
        ALL_TESTS_ARE_CORRECT = false;
    }
    else
    {
        printf("Точка (5 : 1 : 4) не принадлежит кривой\n");
    }


    // ТECT 5: [q]P = E?
    printf("\nТECT 5: [q]P = E?\n");
    struct Point Q;
    mp_int qx, qy, qz;
    mp_init_multi(&qx, &qy, &qz, NULL);
    InitPoint(&Q, &qx, &qy, &qz);

    MontgomeryLadder(&P_base, &ParamSet.q, &Q, &Quadric);
    PrintPoint(&Q);
    PrintPointAffine(&Q, &Quadric);
    if (ArePointsEqual(&Q, &E, &Quadric))
    {
        printf("Точки [q]P и E равны\n");
    }
    else
    {
        printf("Точки [q]P и E не равны\n");
        ALL_TESTS_ARE_CORRECT = false;
    }


    // ТECT 6: [q+1]P = P и [q-1] = -P
    printf("\nТECT 6: [q+1]P = P и [q-1] = -P\n");
    mp_int qminus1, qplus1, buf;
    mp_init_set(&buf, 1);
    mp_init_copy(&qminus1, &ParamSet.q);
    mp_sub(&qminus1, &buf, &qminus1);

    struct Point Q1;
    mp_int qx1, qy1, qz1;
    mp_init_multi(&qx1, &qy1, &qz1, NULL);
    InitPoint(&Q1, &qx1, &qy1, &qz1);

    MontgomeryLadder(&P_base, &qminus1, &Q1, &Quadric);


    mp_init_copy(&qplus1, &ParamSet.q);
    mp_add(&qplus1, &buf, &qplus1);

    struct Point Q2;
    mp_int qx2, qy2, qz2;
    mp_init_multi(&qx2, &qy2, &qz2, NULL);
    InitPoint(&Q2, &qx2, &qy2, &qz2);

    MontgomeryLadder(&P_base, &qplus1, &Q2, &Quadric);

    PrintPointAffine(&Q1, &Quadric);
    // -P = (-PX : PY : PZ)
    struct Point minusP;
    mp_clear(&buf);
    mp_init_copy(&buf, &P_base.X);
    mp_neg_mod(&buf, &buf, &Quadric.p);

    InitPoint(&minusP, &buf, &P_base.Y, &P_base.Z);
    PrintPointAffine(&minusP, &Quadric);
    if (ArePointsEqual(&Q1, &minusP, &Quadric))
    {
        printf("Точки [q - 1]P и -P равны\n");
    }
    else
    {
        printf("Точки [q - 1]P и -P не равны\n");
        ALL_TESTS_ARE_CORRECT = false;
    }

    PrintPointAffine(&Q2, &Quadric);
    PrintPointAffine(&P_base, &Quadric);
    if (ArePointsEqual(&Q2, &P_base, &Quadric))
    {
        printf("Точки [q + 1]P и P равны\n");
    }
    else
    {
        printf("Точки [q + 1]P и P не равны\n");
        ALL_TESTS_ARE_CORRECT = false;
    }


    // ТECT 7: Вычисление [k]P при k = 100
    printf("\nТECT 7: Вычисление [k]P при k = 100\n");
    struct Point kP;
    mp_int kpx, kpy, kpz, k;
    mp_init_set(&k, 100);
    mp_init_multi(&kpx, &kpy, &kpz, NULL);
    InitPoint(&kP, &kpx, &kpy, &kpz);
    MontgomeryLadder(&P_base, &k, &kP, &Quadric);

    PrintPointAffine(&kP, &Quadric);
    if (IsPointOnCurve(&kP, &Quadric))
    {
        printf("Точка [k]P принадлежит кривой\n");
    }
    else
    {
        printf("Точка [k]P не принадлежит кривой\n");
        ALL_TESTS_ARE_CORRECT = false;
    }


    // ТЕСТ 8: Вычисление [k]P для случайного k из диапазона [0, q)
    printf("\nТЕСТ 8: Вычисление [k]P для случайного k из диапазона [0, q)\n");
    mp_prime_rand(&k, 1, 200, NULL);
    mp_to_radix(&k, &answer, 250, NULL, 10);
    printf("k = %s\n", answer);
    MontgomeryLadder(&P_base, &k, &kP, &Quadric);
    PrintPointAffine(&kP, &Quadric);
    if (IsPointOnCurve(&kP, &Quadric))
    {
        printf("Точка [k]P принадлежит кривой\n");
    }
    else
    {
        printf("Точка [k]P не принадлежит кривой\n");
        ALL_TESTS_ARE_CORRECT = false;
    }


    // ТЕСТ 9: [k1]P + [k2]P = [k1 + k2]P?
    printf("\nТЕСТ 9: [k1]P + [k2]P = [k1 + k2]P?\n");
    mp_int k1, k2, k1k2;
    mp_init_multi(&k1, &k2, &k1k2, NULL);

    mp_prime_rand(&k1, 1, 100, NULL);
    mp_to_radix(&k1, &answer, 250, NULL, 10);
    printf("k1 = %s\n", answer);

    mp_prime_rand(&k2, 1, 100, NULL);
    mp_to_radix(&k2, &answer, 250, NULL, 10);
    printf("k2 = %s\n", answer);

    mp_add(&k1, &k2, &k1k2);
    mp_to_radix(&k1k2, &answer, 250, NULL, 10);
    printf("k1 + k2 = %s\n", answer);

    struct Point k1P;
    mp_int k1px, k1py, k1pz;
    mp_init_multi(&k1px, &k1py, &k1pz, NULL);
    InitPoint(&k1P, &k1px, &k1py, &k1pz);

    MontgomeryLadder(&P_base, &k1, &k1P, &Quadric);
    PrintPointAffine(&k1P, &Quadric);
    if (IsPointOnCurve(&k1P, &Quadric))
    {
        printf("Точка [k1]P принадлежит кривой\n");
    }
    else
    {
        printf("Точка [k1]P не принадлежит кривой\n");
        ALL_TESTS_ARE_CORRECT = false;
    }

    struct Point k2P;
    mp_int k2px, k2py, k2pz;
    mp_init_multi(&k2px, &k2py, &k2pz, NULL);
    InitPoint(&k2P, &k2px, &k2py, &k2pz);

    MontgomeryLadder(&P_base, &k2, &k2P, &Quadric);
    PrintPointAffine(&k2P, &Quadric);
    if (IsPointOnCurve(&k2P, &Quadric))
    {
        printf("Точка [k2]P принадлежит кривой\n");
    }
    else
    {
        printf("Точка [k2]P не принадлежит кривой\n");
        ALL_TESTS_ARE_CORRECT = false;
    }

    MontgomeryLadder(&P_base, &k1k2, &kP, &Quadric);
    PrintPointAffine(&kP, &Quadric);
    if (IsPointOnCurve(&kP, &Quadric))
    {
        printf("Точка [k1 + k2]P принадлежит кривой\n");
    }
    else
    {
        printf("Точка [k1 + k2]P не принадлежит кривой\n");
        ALL_TESTS_ARE_CORRECT = false;
    }

    Addition(&k1P, &k2P, &k1P, &Quadric);

    if (ArePointsEqual(&k1P, &kP, &Quadric))
    {
        printf("Точки [k1]P + [k2]P и [k1 + k2]P равны\n");
    }
    else
    {
        printf("Точки [k1]P + [k2]P и [k1 + k2]P не равны\n");
        ALL_TESTS_ARE_CORRECT = false;
    }

    if (ALL_TESTS_ARE_CORRECT)
    {
        printf("\nВСЕ ТЕСТЫ УСПЕШНО ПРОЙДЕНЫ\n");
    }
    else
    {
        printf("\nТЕСТИРОВАНИЕ НЕ ПРОЙДЕНО\n");
    }

    printf("\n");


    // Очистка памяти
    ClearPoint(&P_base);
    ClearPoint(&E);
    ClearPoint(&SumPoint);
    ClearPoint(&P2);
    ClearPoint(&Q);
    ClearPoint(&Q1);
    ClearPoint(&Q2);
    ClearPoint(&minusP);
    ClearPoint(&kP);
    ClearPoint(&k1P);
    ClearPoint(&k2P);

    ClearJacobiQuadric(&Quadric);
    ClearParameters(&ParamSet);

    mp_clear_multi(&p, &a, &b, &q, &x_base, &y_base, &theta, NULL);
    mp_clear_multi(&ex, &ey, &ez, NULL);
    mp_clear_multi(&p2x, &p2y, &p2z, NULL);
    mp_clear_multi(&spx, &spy, &spz, NULL);
    mp_clear_multi(&qx, &qy, &qz, NULL);
    mp_clear_multi(&qx1, &qy1, &qz1, NULL);
    mp_clear_multi(&qx2, &qy2, &qz2, &buf, &qminus1, &qplus1, NULL);
    mp_clear_multi(&kpx, &kpy, &kpz, &k, NULL);
    mp_clear_multi(&k1, &k2, &k1k2, NULL);
    mp_clear_multi(&k1px, &k1py, &k1pz, NULL);
    mp_clear_multi(&k2px, &k2py, &k2pz, NULL);
}

int main()
{
    Test();

    return 0;
}
