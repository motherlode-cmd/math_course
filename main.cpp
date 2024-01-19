#include <iostream>
#include <stdio.h>
#include <vector>
#include "math.h"
using namespace std;
/*_____________________________________________________________________________________________________
----------------------------------------КОНСТАНТЫ ЭКСПЕРИМЕНТА-----------------------------------------
*/
double Ro_w = 1000; //				[ kg/m**3 ] 		плотность воды
double Ro_i = 900; //				[ kg/m**3 ] 		плотность льда
double c_w = 4200; //				[ Дж/(kg * C*) ]	Удельная теплоемкость воды
double c_i = 2100; //				[ Дж/(kg * C*) ]	Удельная теплоемкость льда
double L_i = 340;  //				[ Дж/kg ]			Удельная теплота плавления льда
double k1 = 0.00;  //[ ]					Коэффициент для плавления под водой
double k2 = 0.0000005;  //	    	[ ]					Коэффициент для плавления поверхности над водой
/*
 ----------------------------------------ПЕРЕМЕННЫЕ ЭКСПЕРИМЕНТА----------------------------------------
*/
double m_w = 1000; // 				[ kg ] 				масса воды в начале
double m_i = 100; //					[ kg ] 				масса льда в начале kg
double T_w = 30;//					[ C* ] 				темеатуа воды в начале
double T_a = 26;//					[ C* ] 				темеатуа воздуха в начале
double h = pow(m_i/Ro_i, 1/3);//	[ m ]  				высота кубика
double S_all = h * h * 6; //		[ m**2 ]			Площадь поверхности кубика 
double S_under = 4.6 * S_all / 6;// [ m**2 ] 			Площадь поверхности кусочка льда под водой
double S_ = 1.4 * S_all / 6;//		[ m**2 ]			Площадь поверхности кусочка льда над водой
/*______________________________________________________________________________________________________
*/

// double A = (T_w * m_w * c_w - (m_i * L_i))/Ro_i;
// double B = (m_w + m_i) * c_w / Ro_i;
// double Q = (B/c_w * h - h * h * h * h / 4)/(h * h * h + A);
// double C = T_w/L_i + log(m_w/Ro_i);
// double q =

int right_counter = 0;

double error_lec(double V, double T) {
//  T через V
    double expected = (m_w * c_w * T_w - (m_i - V * Ro_i)*L_i)/(c_w * (m_w + m_i - V * Ro_i));
    return (expected - T)/expected;
}

//Функция, возвращающая площадь грани кубика
double getSquare(double V) { 
	return pow(V, 2/3);
}

typedef struct Data {
	double water_mass; //масса воды в начале
	double ice_mass; //масса льда в начале	
	double k1;
	double k2;
	double T_air; //температура воздуха
} Data;

double error_analit(double V, double T,const Data * data) {
    double M = data->water_mass + data->ice_mass;
    double Vs = data->ice_mass / Ro_i;
    //double expT = T_w * (M - Ro_i * Vs)/(M - Ro_i * V);
    double C = log(T_w * (M - Ro_i * Vs));
    double expT = exp(C) / (M - Ro_i * V);
    return (T - expT)/expT;
    //return (T - C + logc(p - q * V)/q)/T;
}
/*__________________________________________________________________________________
---------Выведенная формула производной для объема льда, отрицательное значение ----
____________________________________________________________________________________
*/

double f_V(double V_current, double T_w_current,const Data * data) 
{
	return -data->k1 * T_w_current * 4.6 * getSquare(V_current) - data->k2 * data->T_air * 1.4 * getSquare(V_current);
}
/*__________________________________________________________________________
----------Выведенная формула прозводной для температуры воды----------------
____________________________________________________________________________
*/

double f_T(double V_current, double T_w_current, const Data * data)
 {	
     double M = data->water_mass + data->ice_mass;
     double Vs = data->ice_mass / Ro_i;
     double C = log(T_w * (M - Ro_i * Vs));
     return f_V(T_w_current, V_current, data) * Ro_i * exp(C)/((M - Ro_i * V_current) * (M - Ro_i * V_current));
     //return f_V(T_w_current, V_current, data) * Ro_i * T_w_current / (M - Ro_i * V_current);
     //return f_V(T_w_current, V_current, data) * (Ro_i * L_i)/
       //                                         (c_w * (data->water_mass + data->ice_mass - V_current * Ro_i));
// 	return (f_V(T_w_current, V_current, data)) * (L_i * Ro_i) /* + Ro_w * T_w_current)*/
// /*	        ___________________________________________________    */
// 			/ (data->water_mass + (data->ice_mass - (V_current * Ro_i) - f_V(T_w_current, V_current, data) * Ro_i));
// /*			           	/						|
//                 масса воды 			 масса льда в начале     текущий объем льда
// */
}
/*_______________________________________________________________________________________________
----------Формула, используемая в пошаговых методах, записывающая значения производных в массив f
_________________________________________________________________________________________________
*/
void f(int n, const double * X, double * f,const Data * data) 
{
	f[0] = f_V(X[0], X[1], data); //	Производная объема
	f[1] = f_T(X[0], X[1], data); //	Производная температуры воды
    //printf("\n v' = %lf T' = %lf\n", f[0], f[1]);
    right_counter++;
}
/*_______________________________________________________________________________________________
----------Метод Эйлера для проверки на корректность работы---------------------------------------
_________________________________________________________________________________________________
*/

void euler_(int n, double * x, void (*f)(int, const double *, double *,const Data *), double h, const Data * data)
{
    double *fx = new double[n];
    f(n,x,fx, data);
    for(int i = 0; i < n; i++) {
        x[i] += fx[i] * h;
    }
    delete [] fx;
}

/*_______________________________________________________________________________________________
----------Таблица Бутчера для метода Ралстона 4-го порядка---------------------------------------
_________________________________________________________________________________________________
*/
double a_ralston[4][4] = {
	{0,		0, 	 0, 	0}, //
	{2.0/5, 0, 	 0, 	0},
	{(-2889 + 1428 * sqrt(5))/1024, (3875 - 1620 * sqrt(5))/1024,0, 0},
	{(-3365 + 2094 * sqrt(5))/6040, (-975 - 3046 * sqrt(5))/2552, (467040+203968*sqrt(5))/240845, 0}
};

double b_ralston[4] = {
	(263 + 24 * sqrt(5))/1812, 
		(125-1000*sqrt(5))/3828,
			(1024*(3346+1623*sqrt(5)))/5924787,
				(30-4*sqrt(5))/123
};

/*_______________________________________________________________________________________________
----------Метод Ралстона 4-го порядка------------------------------------------------------------
_________________________________________________________________________________________________
*/
void ralston(int n, double * x, void (*f)(int, const double *, double *,const Data *), double h, const Data * data)
//             /      |              |                                                   |                |
//            /  переменные(V,T)     |                                                  шаг           Константы эксперимента
{//кол-во переменных              Функция для взятия производной по времени и обьему

    double ** k = new double * [4];
    for(int i = 0; i < 4; i++) {
        k[i] = new double[n];
        double *xt = new double[n];
        for(int pos = 0; pos < n; pos++) {
            xt[pos] = x[pos];
        }
        for(int j = 0; j < 4; j++){
            if(a_ralston[i][j] != 0) {
                for(int pos = 0; pos < n; pos++) {
                    xt[pos] += k[j][pos] * a_ralston[i][j];
                }
            }
        }
        f(n, xt , k[i], data);
        delete [] xt;
    }
    for(int n_ = 0; n_ < n; n_++) {
        for(int i = 0; i < 4; i++) {
            x[n_] += h * b_ralston[i] * k[i][n_];

        }
    }
    for(int i = 0; i < 4; i++) {
        delete [] k[i];
    }
    delete [] k;
}

/*_______________________________________________________________________________________________
----------Неявный метод Гаусса-Лежанда 4-го порядка----------------------------------------------
_________________________________________________________________________________________________
*/
double a_gauss[2][2] = {
	{0.25,     0.25 - sqrt(3)/6}, //
	{0.25 + sqrt(3)/6,     0.25}
};

double b_gauss[2] = {
	0.5, 0.5
};

//x(x0 + h) = x0 + b0*h*k1 + b1 * h * k2;
//k1 = f(x1) = f(x0 + h(a00 * k1 + ao1 * k2))
//k2 = f(x2) = f(x0 + h(a10*k1 + a11*k2)

void gauss_leg_4_ord(int n, double * x, void (*f)(int, const double *, double *,const Data *), double h, const Data * data) {
    double * k1 = new double[n];
    double * k2 = new double[n];
    double * prev_k1 = new double[n];
    double * prev_k2 = new double[n];
    double eps = 0.00001;
    for(int i = 0; i < n; i++) {
        prev_k1[i] = x[i] * 1 + h * (a_gauss[0][0] + a_gauss[0][1]);
        prev_k2[i] = x[i] * 1 + h * (a_gauss[1][0] + a_gauss[1][1]);
    }
    f(n, prev_k1, k1, data);
    f(n, prev_k2, k2, data);
    while(1) {
        for(int i = 0; i < n; i++) {
            prev_k1[i] = k1[i];
            prev_k2[i] = k2[i];
        }
        double * tmp_arg_k1 = new double[n];
        double * tmp_arg_k2 = new double[n];
        for(int i = 0; i < n; i++) {
            tmp_arg_k1[i] = x[i] + h * (k1[i] * a_gauss[0][0] + k2[i] * a_gauss[0][1]);
            tmp_arg_k2[i] = x[i] + h * (k1[i] * a_gauss[1][0] + k2[i] * a_gauss[1][1]);
        }
        f(n, tmp_arg_k1, k1, data);
        f(n, tmp_arg_k2, k2, data);
        delete [] tmp_arg_k1;
        delete [] tmp_arg_k2;
        if(fabs(k1[0] - prev_k1[0]) < eps && fabs(k1[1] - prev_k1[1]) < eps && fabs(k2[0] - prev_k2[0]) < eps && fabs(k2[1] - prev_k2[1]) < eps) {
            break;
        }
    }

    for(int i = 0; i < n; i++) {
        x[i] += b_gauss[0] * h * k1[i] + b_gauss[1] * h * k2[i];
    }
    delete [] prev_k1;
    delete [] prev_k2;
    delete [] k1;
    delete [] k2;
}

/*_______________________________________________________________________________________________
----------Неявный метод Адамса 5-го порядка----------------------------------------------
_________________________________________________________________________________________________
*/
void n_adams5(int n, double* x, double* temp_1, double* temp_2, double* temp_3, double* temp_4, void (*f)(int, const double *, double *,const Data *), double h, const Data * data){
    double * temp_x = new double[n];
    double * temp_5 = new double[n];
    for(int i = 0; i < n; i++) {
        temp_x[i] = x[i] + 55.0/24 * temp_4[i] - 59.0/24 * temp_3[i] + 37.0/24 * temp_2[i] - 3.0/8 * temp_1[i];
    }
    f(n, temp_x, temp_5, data);
    //correct
    for(int i = 0; i < n; i++) {
        temp_x[i] = x[i] + h * (251.0/720 * temp_5[i] + 646.0/720 * temp_4[i] - 264.0/720 * temp_3[i] + 106.0/720 * temp_2[i] - 19.0/720 * temp_1[i]);
    }

    f(n, temp_x, temp_5, data);

    for(int i = 0; i < n; i++) {
        x[i] += h * (251.0/720 * temp_5[i] + 646.0/720 * temp_4[i] - 264.0/720 * temp_3[i] + 106.0/720 * temp_2[i] - 19.0/720 * temp_1[i]);
    }
    for(int i = 0; i < n; i++) {
        temp_1[i] = temp_2[i];
        temp_2[i] = temp_3[i];
        temp_3[i] = temp_4[i];
        temp_4[i] = temp_5[i];
    }
    delete [] temp_x;
    delete [] temp_5;
}

/*_______________________________________________________________________________________________
----------Проверки работ методов(вывод каждого в консоль)----------------------------------------
_________________________________________________________________________________________________
*/
void check_ralston(double * x, double h, double &t, Data & data)
{
    int c = 0;
    while(x[0] > 0 && x[1] > 0)
    {
        //printf("%1lf %1.7lf %2.2lf\n", t, x[0], x[1]);
        //if(c % 50 == 0)
        printf("error %2.4lf %le \n",t,  error_analit(x[0],x[1], &data));
        //printf("%1lf %1.7lf %2.2lf\n", t, x[0], error_lec(x[0],x[1]));

        ralston(2, x, f ,h, &data);
        //euler_(2, x, f ,h, &data);
        t+=h;
        c++;
    }
}

void check_gauss(double * x, double h, double &t, Data & data)
{
    while(x[0] > 0 && x[1] > 0)
    {
        //printf("%1lf %1.7lf %2.6lf\n", t, x[0], x[1]);
        printf("error %2.4lf %le \n",t,  error_analit(x[0],x[1], &data));
        //printf("%1lf %1.7lf %2.2lf\n", t, x[0], error_lec(x[0],x[1]));
        gauss_leg_4_ord(2, x, f, h, &data);
        t+=h;
    }
}

void check_adams(double * x, double h, double &t, Data & data)
{
    int n = 2;
    double * temp_1 = new double[2];
    ralston(2, x, f, h, &data);
    f(n, x, temp_1, &data);
    t+=h;
    //printf("error %2.4lf %le \n",t,  error_analit(x[0],x[1], &data));
    //printf("%1lf %1.5lf %2.2lf\n", t, x[0], x[1]);
    //printf("%1lf %1.7lf %2.2lf\n", t, x[0], error_lec(x[0],x[1]));

    double * temp_2 = new double[2];
    ralston(2, x, f, h, &data);
    f(n, x, temp_2, &data);
    t+=h;
    //printf("error %2.4lf %le \n",t,  error_analit(x[0],x[1], &data));
    //printf("%1lf %1.5lf %2.2lf\n", t, x[0], x[1]);
    //printf("%1lf %1.7lf %2.2lf\n", t, x[0], error_lec(x[0],x[1]));

    double * temp_3 = new double[2];
    ralston(2, x, f, h, &data);
    f(n, x, temp_3, &data);
    t+=h;
    // printf("error %2.4lf %le \n",t,  error_analit(x[0],x[1], &data));
    //printf("%1lf %1.5lf %2.2lf\n", t, x[0], x[1]);
    //printf("%1lf %1.7lf %2.2lf\n", t, x[0], error_lec(x[0],x[1]));

    double * temp_4 = new double[2];
    ralston(2, x, f, h, &data);
    f(n, x, temp_4, &data);
    t+=h;

    while(x[0] > 0 && x[1] > 0)
    {
        //printf("%1lf %1.5lf %2.2lf\n", t, x[0], x[1]);
        printf("error %2.4lf %le \n",t,  error_analit(x[0],x[1], &data));
        t+=h;
        n_adams5(2, x, temp_1, temp_2, temp_3, temp_4, f, h, &data);
    }

    delete [] temp_1;
    delete [] temp_2;
    delete [] temp_3;
    delete [] temp_4;
}

void chaeck_and_compare_all(double * x, double h, double &t, Data & data) {

     //check_adams(x, h, t, data);12224

    //double error_adams =  error_analit(x[0], x[1], &data);


    //check_gauss(x, h, t, data);

    //double error_gauss =  error_analit(x[0], x[1], &data);


    //check_ralston(x, h, t, data);24424

//    double error_ralston =  error_analit(x[0], x[1], &data);


    printf("right counter  : %d\n", right_counter);

}

int main() 
{
	double x[2] = {m_i/Ro_i, T_w};
    double h = 1;
    double t = 0;
    Data data = {m_w, m_i, k1, k2, T_a};
    //printf("%1lf %le %le\n", t, x[0], x[1]);
    //check_adams(x, h, t, data);
    //check_gauss(x, h, t, data);
    //check_ralston(x, h, t, data);
    //printf("A = %1lf B = %le Q = %le , C = %le\n ", A, B, Q, C);
    //printf("error analit %1le\n", error_analit(x[0], x[1], &data));
    //printf("right_counter = %d", right_counter );
   // printf("error LEC %1lf\n ", error_lec(x[0], x[1]));
    chaeck_and_compare_all(x,h,t,data);
}

