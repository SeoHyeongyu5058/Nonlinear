# Non-linear equation

**개요**<br>
* Non-linear equation solver<br>
  * [nonlinearSys2]()<br> 


## nonlinearSys2( )<br>

```c
Matrix nonlinearSys2(Matrix Funcs2, Matrix Jacob2, Matrix Z0, double tol);
```
> Parameter<br>
**Funcs2** : function matrix  <br>
**Jacob2** : Funcs2의 기울기 정보 Jacobian matrix <br>
**z0** : 초기값 정보가 있는 벡터
**tol** : 오차값

<br>
<hr>

1. myMatrix.h 안에 선언되어 있다.<br>
2. 다음과 같은 원리를 기반으로 작성되었다.<br>

* A system of non-linear equations:

$$\bar{X}_{k+1}=\bar{X}_k+H_k$$
$$=\bar{X}_k-(F^{'-1})F(\bar{X}_k)$$
$$=\bar{X}_k-J^{-1}F(\bar{X}_k)$$

* Jacobian Matrix:

$$J=F^{'}({\bar{X}_k})=
\begin{equation}
  \begin{bmatrix}
  dF\over{dx_1} & dF\over{dx_2} & \cdots & dF\over{dx_n}\\
  \end{bmatrix}
\end{equation}
$$


## Example <br>
```c++
int main()
{
	double tol = 0.00001;
	double z0[] = { 30, 100,100 };

	Matrix Z = arr2Mat(z0, 3, 1);


	Matrix F = myFunc2(Z);
	Matrix J = myJacob2(Z);


	nonlinearSys2(F, J, Z, tol);

  double th = Z.at[0][0];
  double dx = Z.at[1][0];
  double dy = Z.at[2][0];

  double t[2] = { dx, dy };
  Matrix _t = arr2Mat(t, 2, 1);
  double R[] = { cos(th * PI / 180), - sin(th * PI / 180), sin(th * PI /180), cos(th * PI/180)};
  Matrix _R = arr2Mat(R, 2, 2);
  double p0[2] = { 0,100 };
  Matrix _p0 = arr2Mat(p0, 2, 1);

  Matrix P0new = createMat(2, 1);
  P0new = addMat(multMat(_R, _p0), _t);
  printMat(P0new, "Matrix P0new");


  freeMat(Z);
  return 0; 
}

Matrix myFunc2(Matrix X)
{
	Matrix F = createMat(X.rows, 1);

	double x0 = 0, y0 = 100, x1 = 0, y1 = -100;
  double x_new0 =50, y_new0 = 186.6025, x_new1 = 150, y_new1 = 13.3975;

	F.at[0][0] = x0 * cos(X.at[0][0] * PI / 180) - y0 * sin(X.at[0][0] * PI / 180) + X.at[1][0] - x_new0;
	F.at[1][0] = x0 * sin(X.at[0][0] * PI / 180) + y0 * cos(X.at[0][0] * PI / 180) + X.at[2][0] - y_new0;
	F.at[2][0] = x1 * cos(X.at[0][0] * PI / 180) - y1 * sin(X.at[0][0] * PI / 180) + X.at[1][0] - y_new1;

	return F;
}

Matrix myJacob2(Matrix X)
{
	Matrix J = createMat(X.rows, X.rows);


	double x0 = 0, y0 = 100, x1 = 0, y1 = -100;

	J.at[0][0] = -y0 * cos(X.at[0][0] * PI / 180) - x0 * sin(X.at[0][0] * PI / 180);
	J.at[0][1] = 1;
	J.at[0][2] = 0;
	J.at[1][0] = x0 * cos(X.at[0][0] * PI / 180) - y0 * sin(X.at[0][0] * PI / 180);
	J.at[1][1] = 0;
	J.at[1][2] = 1;
	J.at[2][0] = -y1 * cos(X.at[0][0] * PI / 180) - x1 * sin(X.at[0][0] * PI / 180);
	J.at[2][1] = 1;
	J.at[2][2] = 0;

	return J;
}
```

## Output <br>
```c
matrix Z =
      17.268264
      31.698750
      75.716589

Matrix P0new =
      50.000000
     186.602540
```
## Warning
>_x.rows = _y.rows<br>
-x.cols=-y.cols


## Error Handling
```c
if (_x.rows != _y.rows )
{
	printf("Error: The sizes of matrix _x and _y are not same\n");

	return zeros(_x.rows, _x.cols);
}
```
<br>
