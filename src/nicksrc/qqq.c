long expmod(long a, long b, long n) 
{ 
 int t ; 
 long ax=1, bx, z, z2 ; 
 t = b % 2 ;  
 if (t==1) ax = a ; 
 bx = b/2; 
 if (bx == 0) return ax % n ; 
 z = expmod(a, bx, n) ; 
 z2 = (z*z) % n ; 
 z2 = (ax*z2) % n ; 
 
 return z2 ; 

}

long
nextprime (long num)
// return nextprime >= num
{
  long x, q;
  int t;

  for (x = num;; ++x) {
    q = expmod(2, x-1, x) ;  
    if (q != 1 ) continue ; 
    t = isprime (x);
    if (t == YES)
      return x;
  }
}

int
isprime (long num)
// naive algorithm.  Implement Pollard rho at some time
{
  int top, x, t;

  if (num < 2)
    return NO;
  if (num == 2)
    return YES;
  top = nnint (sqrt (num));

  t = num % 2 ; 
  if (t==0) return NO ;  

  for (x = 3; x <= top; x += 2) {
    t = num % x;
    if (t == 0)
      return NO;
  }

  return YES;

}


