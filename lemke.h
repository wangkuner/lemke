#ifndef _LEMKE_H
#define _LEMKE_H

enum class Base_Type
{
    zero = 0,
    one = 1
};

typedef double _real_lemke;
typedef int _int_lemke;

struct Base_Variable
{
    Base_Type type;
    _int_lemke id;
};

int Lemke(const _int_lemke n, const _real_lemke *M, const _real_lemke *q,
          const _int_lemke max_iter, _int_lemke &num_iter,
          _real_lemke *x, _real_lemke *z, _real_lemke &error);

#endif //_LEMKE_H
