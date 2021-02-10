#ifndef NULL
#define NULL ((void*)0)
#endif


#ifndef CURVE_H
#define CURVE_H
#include "tommath.h"

// Набор параметров из стандарта
#define p_str   "115792089237316195423570985008687907853269984665640564039457584007913129639319"
#define q_str   "28948022309329048855892746252171976963338560298092253442512153408785530358887"

#define a_str   "87789765485885808793369751294406841171614589925193456909855962166505018127157"
#define b_str   "18713751737015403763890503457318596560459867796169830279162511461744901002515"

#define x_base_str  "65987350182584560790308640619586834712105545126269759365406768962453298326056"
#define y_base_str  "22855189202984962870421402504110399293152235382908105741749987405721320435292"

// Значение, вычисленное при помощи Wolfram Mathematica на основе параметров из стандарта
#define theta_str   "454069018412434321972378083527459607666454479745512801572100703902391945898"


// Точка в проективных координатах (X : Y : Z)
struct Point
{
    mp_int X;
    mp_int Y;
    mp_int Z;
};

// Параметры в виде больших чисел
struct Parameters
{
    mp_int p;
    mp_int q;

    mp_int a;
    mp_int b;

    mp_int x_base;
    mp_int y_base;

    mp_int theta;
};

// Параметры квадрики Якоби
struct JacobiQuadric
{
    mp_int p;
    mp_int e;
    mp_int d;

    mp_int X;
    mp_int Y;
    mp_int Z;
};


/*
 * Функции инициализации основных структур
 */
// Инициализация точки
void InitPoint(struct Point * P, mp_int * x, mp_int * y, mp_int * z);

// Инициализация структуры с параметрами
void InitParameters(struct Parameters * Param, mp_int * p, mp_int * q, mp_int * a, mp_int * b, mp_int * x_base, mp_int * y_base, mp_int * theta);

// Инициализация структуры с параметрами кривой
void InitJacobiQuadric(struct JacobiQuadric * JQ, struct Parameters * Param);


/*
 * Функции освобождения памяти, использующейся для хранения основных структур
 */
// Очистка памяти, использовавшейся для хранения точки
void ClearPoint(struct Point * P);

// Очистка памяти, использовавшейся для хранения параметров
void ClearParameters(struct Parameters * Param);

// Очистка памяти, использовавшейся для хранения параметров кривой Якоби
void ClearJacobiQuadric(struct JacobiQuadric * JQ);


/*
 * Функции для вывода координат точки на экран
 */
// Вывести значения проективных координат точки на экран
void PrintPoint(struct Point * P);

// Вывести значения афинных координат точки на экран
void PrintPointAffine(struct Point * P, struct JacobiQuadric * JQ);


/*
 * Основные функции
 */
// Проверить находится ли данная точка на данной кривой
int IsPointOnCurve(struct Point * P, struct JacobiQuadric * JQ);

// Сложение двух точек P1 + P2, результат записывается в третью точку P3
void Addition(struct Point * P1, struct Point * P2, struct Point * P3, struct JacobiQuadric * JQ);

// Проверка равенства двух точек на кривой
int ArePointsEqual(struct Point * P1, struct Point * P2, struct JacobiQuadric * JQ);

// Реализация алгоритма "лесенка Монтгомери"
void MontgomeryLadder(struct Point * P, mp_int * k, struct Point * Q, struct JacobiQuadric * JQ);

#endif
