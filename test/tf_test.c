#include <dzco/dz_tf.h>

#define TOL ( 1.0e-7 ) /* tolerance to compare values */

bool is_equal(double v1, double v2)
{
  bool ret;
  ret = ( zIsTol(v1,TOL) && zIsTol(v2,TOL) ) ||
        zIsTol( v1/v2 - 1.0, TOL ) ? true : false;
  return ret;
}

bool complex_is_equal(zComplex *c1, zComplex *c2)
{
  return is_equal( c1->re, c2->re ) && is_equal( c1->im, c2->im ) ? true : false;
}

bool cval_set_is_included(zCVec v, zComplex *val)
{
  int i;
  bool ret = false;

  for( i=0; i<zCVecSizeNC(v); i++ ){
    if( complex_is_equal( val, zCVecElemNC(v,i) ) ) ret = true;
  }
  return ret;
}

bool cval_set_equal(zCVec v1, zCVec v2)
{
  int i;
  bool ret = true;

  for( i=0; i<zCVecSizeNC(v2); i++ )
    if( !cval_set_is_included( v1, zCVecElemNC(v2,i) ) ) ret = false;
  for( i=0; i<zCVecSizeNC(v1); i++ )
    if( !cval_set_is_included( v2, zCVecElemNC(v1,i) ) ) ret = false;
  return ret;
}

#define NUM_ZEROS 2
#define NUM_POLES 5

void assert_tf_zeropole(void)
{
  zCVec zero_src, pole_src;
  dzTF tf;

  zero_src = zCVecAlloc( NUM_ZEROS );
  pole_src = zCVecAlloc( NUM_POLES );

  zComplexCreate( zCVecElemNC(zero_src,0), zRandF(-10,0), 0 );
  zComplexCreate( zCVecElemNC(zero_src,1), zRandF(0,10), 0 );
  zComplexCreate( zCVecElemNC(pole_src,0), zRandF(-10,0), 0 );
  zComplexCreate( zCVecElemNC(pole_src,1), zRandF(-10,0), zRandF(-10,10) );
  zComplexConj( zCVecElemNC(pole_src,1), zCVecElemNC(pole_src,2) );
  zComplexCreate( zCVecElemNC(pole_src,3), zRandF(-10,0), zRandF(-10,10) );
  zComplexConj( zCVecElemNC(pole_src,3), zCVecElemNC(pole_src,4) );

  dzTFCreateZeroPole( &tf, zero_src, pole_src, zRandF(0.1,5) );
  dzTFZeroPole( &tf );

  zAssert( dzTFCreateZeroPole + dzTFZeroPole,
    cval_set_equal( zero_src, dzTFZero(&tf) ) &&
    cval_set_equal( pole_src, dzTFPole(&tf) ) );

  dzTFDestroy( &tf );
  zCVecFree( zero_src );
  zCVecFree( pole_src );
}

void assert_tf_connect(void)
{
  dzTF tf1, tf2, tf;

  dzTFAlloc( &tf1, 2, 3 );
  dzTFAlloc( &tf2, 2, 3 );
  dzTFSetNumList( &tf1, zRandF(-10,10), zRandF(-10,10), zRandF(-10,10) );
  dzTFSetDenList( &tf1, zRandF(-10,10), zRandF(-10,10), zRandF(-10,10), zRandF(-10,10) );
  dzTFSetNumList( &tf2, zRandF(-10,10), zRandF(-10,10), zRandF(-10,10) );
  dzTFSetDenList( &tf2, zRandF(-10,10), zRandF(-10,10), zRandF(-10,10), zRandF(-10,10) );
  dzTFNum(&tf) = zPexMul( dzTFNum(&tf1), dzTFNum(&tf2) );
  dzTFDen(&tf) = zPexMul( dzTFDen(&tf1), dzTFDen(&tf2) );
  dzTFConnect( &tf1, &tf2 );
  zAssert( dzTFConnect,
    zPexEqual(dzTFNum(&tf),dzTFNum(&tf1),zTOL) &&
    zPexEqual(dzTFDen(&tf),dzTFDen(&tf1),zTOL) );
}

/* zeros and poles of a transfer function with a double pole, a triple zero, and a pair of complex poles,
 * namely, (s+2)^3 / ( (s+1)^2 (s^2+2s+5) ), classified into real and complex values. */
void assert_tf_zeropole_reim(void)
{
  dzTF tf;
  zVec zero1 = NULL, pole1 = NULL;
  zCVec zero2 = NULL, pole2 = NULL;
  int i;
  bool result = true;

  dzTFInit( &tf );
  dzTFAlloc( &tf, 3, 4 );
  dzTFSetNumList( &tf, 8.0, 12.0, 6.0, 1.0 );
  dzTFSetDenList( &tf, 5.0, 12.0, 10.0, 4.0, 1.0 );
  if( !dzTFZeroPoleReImDefault( &tf, &zero1, &zero2, &pole1, &pole2 ) ||
      !zero1 || zVecSizeNC(zero1) != 3 || zero2 ||
      !pole1 || zVecSizeNC(pole1) != 2 || !pole2 || zCVecSizeNC(pole2) != 2 ) result = false;
  for( i=0; result && i<3; i++ )
    if( !zIsTol( zVecElemNC(zero1,i) + 2, 1.0e-4 ) ) result = false;
  for( i=0; result && i<2; i++ )
    if( !zIsTol( zVecElemNC(pole1,i) + 1, 1.0e-6 ) ) result = false;
  zVecFree( zero1 );
  zCVecFree( zero2 );
  zVecFree( pole1 );
  zCVecFree( pole2 );
  dzTFDestroy( &tf );
  zAssert( dzTFZeroPoleReImDefault, result );
}

int main(int argc, char *argv[])
{
  zRandInit();
  assert_tf_zeropole();
  assert_tf_connect();
  assert_tf_zeropole_reim();
  return 0;
}
