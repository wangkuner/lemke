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

using std::memcpy;
using std::memset;

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
          _real_lemke *x, _real_lemke *z, _real_lemke &error)
{
    _int_lemke i = 0;
    _int_lemke j = 0;
    num_iter = 0;
    error = 0.0;
    memset(x, 0, sizeof(x[0]) * n);
    memcpy(z, q, sizeof(x[0]) * n);
    while (i < n && q[i] > 0.0)
    {
        i++;
    }
    if (i >= n)
    {
        return 0;
    }
    _int_lemke s = 0;
    _int_lemke r = 0;
    _int_lemke next_enter = 0;
    const _real_lemke eps = 1e-16;
    _int_lemke col = 2 * n + 2;
    _real_lemke *table = new _real_lemke[n * col];
    _int_lemke *leave_entry = new _int_lemke[n];
    _real_lemke min = 0.0;

    memset(table, 0, sizeof(_real_lemke) * n * col);
    for (i = 0; i < n; i++)
    {
        table[i * col + i] = 1.0;
        for (j = 0; j < n; j++)
        {
            table[i * col + n + j] = -M[i * n + j];
        }
        table[i * col + n + n] = -1.0;
        table[i * col + n + n + 1] = q[i];
        // z index [0,n-1], x index [n,2n-1]. leave_entry[i] is z[i]
        leave_entry[i] = i;
    }
    // column 2*n
    s = 2 * n;
    // line r
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
    if (leave_entry[r] < n)
    {
        next_enter = leave_entry[r] + n; // next enter entry is x[leave_entry[r]]
    }
    else
    {
        next_enter = leave_entry[r] - n; // next enter entry is z[leave_entry[r]-n]
    }
    leave_entry[r] = s;
    if (_Pivot(table, n, col, r, s, eps))
    {
        if (table)
        {
            delete[] table;
        }
        table = nullptr;
        if (leave_entry)
        {
            delete[] leave_entry;
        }
        leave_entry = nullptr;
        return 1;
    }
    for (num_iter = 0; num_iter < max_iter; ++num_iter)
    {
        // find column s
        s = next_enter;
        // find line r
        i = 0;
        while (i < n && table[i * col + s] < eps)
        {
            ++i;
        }
        // all elements are non-positive
        if (i >= n || table[i * col + s] < eps)
        {
            break;
        }
        min = table[i * col + n + n + 1] / table[i * col + s];
        r = i;
        for (i = i + 1; i < n; i++)
        {
            if (table[i * col + s] > eps)
            {
                if (min > table[i * col + n + n + 1] / table[i * col + s])
                {
                    min = table[i * col + n + n + 1] / table[i * col + s];
                    r = i;
                }
            }
        }
        //the leave base is x0
        if (leave_entry[r] == n + n)
        {
            leave_entry[r] = s;
            if (_Pivot(table, n, col, r, s, eps))
            {
                if (table)
                {
                    delete[] table;
                }
                table = nullptr;
                if (leave_entry)
                {
                    delete[] leave_entry;
                }
                leave_entry = nullptr;
                return 1;
            }
            break;
        }
        else //leave z_r. enter y_s.
        {
            if (leave_entry[r] < n)
            {
                next_enter = leave_entry[r] + n; // next enter entry is x[leave_entry[r]]
            }
            else
            {
                next_enter = leave_entry[r] - n; // next enter entry is z[leave_entry[r]-n]
            }
            leave_entry[r] = s;
            if (_Pivot(table, n, col, r, s, eps))
            {
                if (table)
                {
                    delete[] table;
                }
                table = nullptr;
                if (leave_entry)
                {
                    delete[] leave_entry;
                }
                leave_entry = nullptr;
                return 1;
            }
        }
    }
    memset(z, 0, sizeof(x[0]) * n);
    for (i = 0; i < n; i++)
    {
        if (leave_entry[i] < n)
        {
            z[leave_entry[i]] = table[i * col + n + n + 1];
        }
        else if (leave_entry[i] < n + n)
        {
            x[leave_entry[i] - n] = table[i * col + n + n + 1];
        }
    }
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
    if (leave_entry)
    {
        delete[] leave_entry;
    }
    leave_entry = nullptr;

    return 0;
}