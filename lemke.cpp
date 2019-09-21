/*------------------------------------------------------------------------------
*
* 	@Author:Wang Kun
* 	@Data:2019-09-21
*   @Description: Solve the linear complementrary problem.
*                 z = Mx + q, x>=0, z>=0, <x,z>=0.
*                 error=|| <(Mx+q),x> ||
*                 return the result x, error and iteration number.
*   @References: Yassine, Adnan. Comparative study between lemke's 
*                method and the interior point method for the monotone 
*                linear complementary problem. 2008.01.
------------------------------------------------------------------------------*/

#include "lemke.h"

#include <cstring>

static int _Pivot(_real_lemke *M, const _int_lemke m, const _int_lemke n,
                  const _int_lemke i, const _int_lemke j, const _real_lemke eps)
{
    if (M[i * n + j] < eps && M[i * n + j] > -eps)
    {
        return 1;
    }
    _int_lemke k = 0;
    _int_lemke s = 0;
    _real_lemke pivot = M[i * n + j];
    for (k = 0; k < n; k++)
    {
        M[i * n + k] /= pivot;
    }
    for (s = 0; s < m; s++)
    {
        if (s != i)
        {
            pivot = M[s * n + j];
            for (k = 0; k < n; k++)
            {
                M[s * n + k] -= pivot * M[i * n + k];
            }
        }
    }
    return 0;
}

int Lemke(const _int_lemke n, const _real_lemke *M, const _real_lemke *q,
          const _int_lemke max_iter, _int_lemke &num_iter,
          _real_lemke *x, _real_lemke &error)
{
    _int_lemke i = 0;
    while (i < n && q[i] > 0.0)
    {
        i++;
    }
    if (i >= n)
    {
        return 0;
    }
    bool all_negative = false;
    _int_lemke s = 0;
    _int_lemke r = 0;
    _int_lemke step = 0;
    const _real_lemke eps = 1e-16;
    _real_lemke *table = new _real_lemke[n * (n + 2)];
    _real_lemke min = 0.0;
    Base_Variable(*base)[2] = new Base_Variable[n + 1][2];
    Base_Variable(*xz)[2] = new Base_Variable[n + 1][2];
    Base_Variable record;
    Base_Variable temp;
    for (i = 0; i < n; i++)
    {
        table[i * (n + 2)] = -1.0;
        for (s = 0; s < n; s++)
        {
            table[i * (n + 2) + s + 1] = -M[i * n + s];
        }
        table[i * (n + 2) + n + 1] = q[i];
        base[i][0].type = Base_Type::zero;
        base[i][0].id = i;
        xz[i][0].type = Base_Type::zero;
        xz[i][0].id = i;
        base[i][1].type = Base_Type::one;
        base[i][1].id = i;
        xz[i][1].type = Base_Type::one;
        xz[i][1].id = i;
    }
    base[n][0].type = Base_Type::zero;
    base[n][0].id = n;
    xz[n][0].type = Base_Type::zero;
    xz[n][0].id = n;
    min = q[0];
    r = 0;
    for (i = 1; i < n; i++)
    {
        if (min > q[i])
        {
            min = q[i];
            r = i;
        }
    }
    record.id = 0; //record y_s=x_0
    record.type = Base_Type::zero;
    //leave z_r, enter x_0
    temp = base[r][1];
    base[r][1] = record;
    base[xz[record.id][static_cast<_int_lemke>(record.type)].id][0] = temp;
    //record the new position of x and z
    xz[temp.id][static_cast<_int_lemke>(temp.type)].type = Base_Type::zero;
    xz[temp.id][static_cast<_int_lemke>(temp.type)].id = xz[record.id][static_cast<_int_lemke>(record.type)].id;
    xz[record.id][static_cast<_int_lemke>(record.type)].type = Base_Type::one;
    xz[record.id][static_cast<_int_lemke>(record.type)].id = r;
    record.id = temp.id + 1; //record y_s=x_r
    record.type = static_cast<Base_Type>(1 - static_cast<_int_lemke>(temp.type));
    if (_Pivot(table, n, n + 2, r, 0, eps))
    {
        return 1;
    }
    for (step = 0; step < max_iter; ++step)
    {
        //find column s
        s = xz[record.id][static_cast<_int_lemke>(record.type)].id;
        //find line r
        i = 0;
        while (i < n && table[i * (n + 2) + s] < 0.0)
        {
            ++i;
        }
        if (i >= n || table[i * (n + 2) + s] < eps)
        {
            return 1;
        }
        min = table[i * (n + 2) + n + 1] / table[i * (n + 2) + s];
        r = i;
        for (i = i + 1; i < n; i++)
        {
            if (table[i * (n + 2) + s] > eps)
            {
                if (min > table[(i + 1) * (n + 2) - 1] / table[i * (n + 2) + s])
                {
                    min = table[(i + 1) * (n + 2) - 1] / table[i * (n + 2) + s];
                    r = i;
                }
            }
        }
        //the enter base is x0
        if (base[r][1].type == Base_Type::zero && base[r][1].id == 0)
        {
            if (_Pivot(table, n, n + 2, r, s, eps))
            {
                return 1;
            }
            base[r][1] = record;
            break;
        }
        else //leave z_r. enter y_s.
        {
            temp = base[r][1];
            base[r][1] = record; //leave z_r, enter x_0
            base[xz[record.id][static_cast<_int_lemke>(record.type)].id][0] = temp;
            xz[temp.id][static_cast<_int_lemke>(temp.type)].type = Base_Type::zero;
            xz[temp.id][static_cast<_int_lemke>(temp.type)].id = xz[record.id][static_cast<_int_lemke>(record.type)].id;
            xz[record.id][static_cast<_int_lemke>(record.type)].type = Base_Type::one;
            xz[record.id][static_cast<_int_lemke>(record.type)].id = r;
            record.id = temp.id + 1; //record y_s=x_s
            record.type = static_cast<Base_Type>(1 - static_cast<_int_lemke>(temp.type));
            if (_Pivot(table, n, n + 2, r, s, eps))
            {
                return 1;
            }
        }
    }
    memset(x, 0, sizeof(x[0]) * n);
    for (i = 0; i < n; i++)
    {
        if (base[i][1].type == Base_Type::zero)
        {
            x[base[i][1].id - 1] = table[(i + 1) * (n + 2) - 1];
        }
    }
    num_iter = step;
    error = 0.0;
    for (i = 0; i < n; i++)
    {
        error += x[i] * q[i];
        for (s = 0; s < n; s++)
        {
            error += x[i] * M[i * n + s] * x[s];
        }
    }

    if (table)
    {
        delete[] table;
    }
    table = nullptr;
    if (xz)
    {
        delete[] xz;
    }
    xz = nullptr;
    if (base)
    {
        delete[] base;
    }
    base = nullptr;
    return 0;
}