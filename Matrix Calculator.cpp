#include<stdio.h>
#include<math.h>
#include<string.h>
#include<stdlib.h>

const int M=20,N=20;
int choice,n,m,s;
double mtx1[M][N],mtx2[N][M],mtx3[M][N];

void MATRIX_INPUT(double (*a)[M]){
	printf("请输入行数：");
	scanf("%d",&m);
	printf("请输入列数：");
	scanf("%d",&n);
	double t;
	for(int i=0;i<m;i++){
		printf("请输入第%d行：",i+1);
		for (int j = 0; j < n; j++){
			scanf("%lf",&t);
			*(*(a+i)+j)=t;
		}
	}
}

void MATRIX_OUTPUT(double (*a)[M]){
	for(int i=0;i<m;i++){
		printf("\n");
		for(int j=0;j<n;j++){
			if(*(*(a+i)+j)==0) *(*(a+i)+j)=0;
			if(*(*(a+i)+j)==(int)*(*(a+i)+j))
				printf("%-5.0f",*(*(a+i)+j));
			else printf("%-5.3f",*(*(a+i)+j));
		}
	}
}

void MATRIX_ROW_ECHELON(double (*a)[M],double (*at)[M],int m,int n){
	int tmp,k=0;
	double ratio;
	for(int i=0;i<m;i++)
		for(int j=0;j<n;j++)
			*(*(at+i)+j)=*(*(a+i)+j);
	for(int i=0;i<m;i++){
		while(k<n){
			for(int q=i;q<m;q++)
				if(*(*(at+q)+k)!=0){
					for(int r=0;r<n;r++){
						int t=*(*(at+i)+r);
						*(*(at+i)+r)=*(*(at+q)+r);
						*(*(at+q)+r)=t;
					}
					goto out;
				}
			k++;
		}
		out:
		if(k==n) return;
		for(int j=i+1;j<m;j++){
			ratio=*(*(at+j)+k)/ *(*(at+i)+k);
			int p=k;
			while(p<n){
				*(*(at+j)+p)-=ratio* *(*(at+i)+p);
				p++;
			}
		}
	}
}

void MATRIX_ROW_SIMPLIFY(double (*a)[M],double (*at)[M],int m,int n){
	MATRIX_ROW_ECHELON(a,at,m,n);
	double ratio1,ratio2;
	for(int i=m-1;i>=0;i--){
		int index1=-1;
		for(int j=0;j<n;j++){
			if(*(*(at+i)+j)!=0 && index1==-1){
				index1=j;
				ratio1=*(*(at+i)+index1);
			}
			if(index1>=0) *(*(at+i)+j)/=ratio1;
		}
		if(index1>=0) for(int k=i-1;k>=0;k--){
			ratio2=*(*(at+k)+index1)/ *(*(at+i)+index1);
			for(int p=index1;p<n;p++)
				*(*(at+k)+p)-=ratio2* *(*(at+i)+p);
		}
	}
}

void MATRIX_TRANSPOSITION(double (*a)[N],double (*at)[M],int m,int n){
	for(int i=0;i<m;i++)
		for(int j=0;j<n;j++)
			*(*(at+i)+j)=*(*(a+j)+i);
}

double MATRIX_DETERMINANT(double **a,int x,int y){
	if(y==1) return a[0][0];
	double det=0;
	int i,j,k,offset = 0;
	double **b=(double**)malloc(sizeof(double)*x);
	for(int i=0;i<x;i++)
		b[i]=(double*)malloc(sizeof(double)*y);
	for(int k=0;k<y;i++){
		for(int i=0;i<x-1;i++)
			for (int j=0;j<y;j++){
				if(j<k) offset=0;
				else if(j>=k) offset=1;
				b[i][j]=*(*(a+i+1)+j+offset);
			}
		printf("x-1=%d y-1=%d\n",x-1,y-1);
		det=det+*(*a+k)*pow(-1,k)*MATRIX_DETERMINANT(b,x-1,y-1);
		printf("det=%d\n",det);
	}
	return det;
}

void MATRIX_ADD(double (*a1)[M],double (*a2)[M],double (*at)[M],int m,int n){
	for(int i=0;i<m;i++)
		for(int j=0;j<n;j++)
			*(*(at+i)+j)=*(*(a1+i)+j)+*(*(a2+i)+j);
}

void MATRIX_MULTIPLY_NUM(int number,double (*a)[M],double (*at)[M],int m,int n){
	int i,j;
	for(i=0;i<m;i++)
		for(j=0;j<n;j++)
			*(*(at+i)+j)=*(*(a+i)+j)*number;
}

void MATRIX_MULTIPLY_MATRIX(double (*a1)[M],double (*a2)[M],double (*at)[M],int m,int s,int n){
	int i,j,k;
	for(i=0;i<m;i++)
		for(j=0;j<n;j++){
			*(*(at+i)+j)=0;
			for(k=0; k < s; k++)
				*(*(at+i)+j)+=*(*(a1+i)+k) * *(*(a2+k)+j);
		}
}

int MATRIX_RANK(double (*a)[M],double (*at)[M],int m,int n){
	int i,j,rank=0;
	MATRIX_ROW_ECHELON(a,at,m,n);
	for(i=0;i<m;i++)
		for(j=0;j<n;j++)
			if(*(*(at+i)+j)!=0){
				rank++;
				break;
			}
	return rank;
}

int main(){
	int x,y;
	printf("Matrix Calculator\n1 求行阶梯形矩阵\n\2 求行最简形矩阵\n3 求行列式的值\n4 矩阵加法\n5 矩阵数乘\n6 矩阵乘法\n7 矩阵转置\n8 求伴随矩阵\n9 矩阵求逆\n10 矩阵求秩\n11 解线性方程组\n");
	scanf("%d",&choice);
	switch(choice){
		case 1:
			MATRIX_INPUT(mtx1);
			MATRIX_ROW_ECHELON(mtx1,mtx2,m,n);
			MATRIX_OUTPUT(mtx2);
			break;
		case 2:
			MATRIX_INPUT(mtx1);
			MATRIX_ROW_SIMPLIFY(mtx1,mtx2,m,n);
			MATRIX_OUTPUT(mtx2);
			break;
		case 3:
			MATRIX_INPUT(mtx1);
			double det;
			x=m;y=n;
			//det = MATRIX_DETERMINANT(mtx1, x, y);
			printf("矩阵的行列式是%f\n",det);
			break;
		case 4:
			printf("请输入第一个矩阵：\n");
			MATRIX_INPUT(mtx1);
			printf("请输入第二个矩阵：\n");
			MATRIX_INPUT(mtx2);
			MATRIX_ADD(mtx1,mtx2,mtx3,m,n);
			MATRIX_OUTPUT(mtx3);
			break;
		case 6:
			printf("请输入第一个矩阵：\n请输入行数：");
			scanf("%d",&m);
			printf("请输入列数：");
			scanf("%d",&s);
			for(int i=0;i<m;i++){
				printf("请输入第%d行：",i+1);
				for(int j=0;j<s;j++)
					scanf("%lf",&mtx1[i][j]);
			}
			printf("请输入第二个矩阵：\n请输入行数：");
			scanf("%d",&s);
			printf("请输入列数：");
			scanf("%d",&n);
			for(int i=0;i<s;i++){
				printf("请输入第%d行：",i+1);
				for(int j=0;j<n;j++)
					scanf("%lf",&mtx2[i][j]);
			}
			MATRIX_MULTIPLY_MATRIX(mtx1,mtx2,mtx3,m,s,n);
			MATRIX_OUTPUT(mtx3);
			break;
		case 7:
			MATRIX_INPUT(mtx1);
			MATRIX_TRANSPOSITION(mtx1,mtx2,m,n);
			MATRIX_OUTPUT(mtx2);
			break;
		case 10:
			MATRIX_INPUT(mtx1);
			int rank;
			rank=MATRIX_RANK(mtx1,mtx2,m,n);
			printf("矩阵的秩是%d\n",rank);
			break;
		default:
			break;
	}
	return 0;
}




















