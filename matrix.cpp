matrix::matrix()
{
    m=n=0;
}

matrix::~matrix()
{
    for(int i=0;i<m;i++) delete[] e[i];
    if(n!=0) {delete[] e;}
    //cout<<"A "<<m<<"x"<<n<<" matrix has been destroyed."<<endl;
}

void matrix::init(int a)
{
    m=a;
    n=a;
    e=new double*[m];
    for(int i=0;i<m;i++) e[i]=new double[n];
    for(int i=0;i<m;i++) {for(int j=0;j<n;j++) {e[i][j]=0.;}}
    //cout<<"A "<<m<<"x"<<n<<" matrix has been created."<<endl;
}

void matrix::init(int a,int b)
{
    m=a;
    n=b;
    e=new double*[m];
    for(int i=0;i<m;i++) e[i]=new double[n];
    for(int i=0;i<m;i++) {for(int j=0;j<n;j++) {e[i][j]=0.;}}
    //cout<<"A "<<m<<"x"<<n<<" matrix has been created."<<endl;
}

void matrix::init(string fn)
{
    double x;
    string s;
    ifstream in;
    in.open(fn.c_str());
    while(in.peek()!='\n') {in>>x;n++;} in.seekg(0,ios::beg);
    while(!in.eof()) {getline(in,s);m++;} m--; in.clear(); in.seekg(0,ios::beg);
    e=new double*[m];
    for(int i=0;i<m;i++) e[i]=new double[n];
    //cout<<"A "<<m<<"x"<<n<<" matrix has been created."<<endl;
    for(int i=0;i<m;i++) {for(int j=0;j<n;j++) {in>>e[i][j];}}
    in.close();
}

void matrix::ident()
{
    for(int i=0;i<m;i++) {for(int j=0;j<n;j++) {if(i==j){e[i][j]=1.;}else{e[i][j]=0.;}}}
}

istream& operator>>(istream &min, matrix &a)
{
    double x;
    string s;
    while(min.peek()!='\n') {min>>x;a.n++;} min.seekg(0,ios::beg);
    while(!min.eof()) {getline(min,s);a.m++;} a.m--; min.clear(); min.seekg(0,ios::beg);
    a.e=new double*[a.m];
    for(int i=0;i<a.m;i++) a.e[i]=new double[a.n];
    //cout<<"A "<<a.m<<"x"<<a.n<<" matrix has been created."<<endl;
    for(int i=0;i<a.m;i++) {for(int j=0;j<a.n;j++) {min>>a.e[i][j];}}
    return min;
}

ostream& operator<<(ostream &mout, matrix &a)
{
    for(int i=0;i<a.m;i++) {for(int j=0;j<a.n;j++) {mout<<a.e[i][j]<<" ";}mout<<endl;}
    return mout;
}

matrix& matrix::operator=(const double &x)
{
    for(int i=0;i<m;i++) {for(int j=0;j<n;j++) e[i][j]=x;} return *this;
}

matrix& matrix::operator=(const matrix &a)
{
    if((m==0)&&(n==0)) {m=a.m;n=a.n;init(m,n);}
    if((m!=a.m)||(n!=a.n)) {cout<<"Matrix Assingment Failure : Dimensional Inequality."<<endl;exit(1);}
    for(int i=0;i<m;i++) {for(int j=0;j<n;j++) {e[i][j]=a.e[i][j];}} return (*this);
}

matrix matrix::operator+(matrix &a)
{
    if((m!=a.m)||(n!=a.n)) {cout<<"Matrix Addition Failure : Dimensional Inequality."<<endl;exit(1);}
    matrix s;
    s.init(m,n);
    for(int i=0;i<m;i++) {for(int j=0;j<n;j++) s.e[i][j]=e[i][j]+a.e[i][j];}
    return s;
}

matrix matrix::operator*(double x)
{
    matrix s;
    s.init(m,n);
    for(int i=0;i<m;i++) {for(int j=0;j<n;j++) s.e[i][j]=x*e[i][j];}
    return s;
}

matrix matrix::operator*(matrix &a)
{
    if(n!=a.m) {cout<<"Matrix Product Failure : Column-Row Inequality."<<endl;exit(1);}
    matrix s;
    s.init(n,a.n);
    for(int i=0;i<m;i++) {for(int j=0;j<a.n;j++) {for(int k=0;k<n;k++) {s.e[i][j]+=e[i][k]*a.e[k][j];}}}
    return s;
}

matrix operator*(double x,matrix &a)
{
    matrix s;
    s=a*x;
    return s;
}

matrix matrix::operator-(matrix &a)
{
    if((m!=a.m)||(n!=a.n)) {cout<<"Matrix Subtraction Failure : Dimensional Inequality."<<endl;exit(1);}
    matrix s;
    s=(a*(-1.))+(*this);
    return s;
}

matrix matrix::operator/(double x)
{
    matrix s;
    s=(*this)*(1./x);
    return s;
}

matrix& matrix::operator+=(matrix &a) {return (*this=*this+a);}
matrix& matrix::operator-=(matrix &a) {return (*this=*this-a);}
matrix& matrix::operator*=(double x) {return (*this=*this*x);}
matrix& matrix::operator*=(matrix &a) {return (*this=*this*a);}
matrix& matrix::operator/=(double x) {return (*this=*this/x);}

bool matrix::operator==(matrix &a)
{
    if((m!=a.m)||(n!=a.n)) {cout<<"Matrix Equality Comparison Failure : Dimensional Inequality."<<endl;exit(1);}
    bool c=true;
    for(int i=0;i<m;i++) {for(int j=0;j<n;j++) {if(e[i][j]!=a.e[i][j]) {c=false;break;}}}
    return c;
}

bool matrix::operator!=(matrix &a)
{
    if((m!=a.m)||(n!=a.n)) {cout<<"Matrix InEquality Comparison Failure : Dimensional Inequality."<<endl;exit(1);}
    bool c=true;
    if(*this==a) {c=false;}
    return c;
}

bool matrix::operator==(double x)
{
    bool c=true;
    for(int i=0;i<m;i++) {for(int j=0;j<n;j++) {if(e[i][j]!=x) {c=false;break;}}}
    return c;
}

bool matrix::operator!=(double x)
{
    bool c=true;
    if(*this==x) {c=false;}
    return c;
}

double matrix::Norm(char L)
{
    double s=0.;
    matrix x;
    x.init(1,m);
    switch(L)
    {
	  case '1' :
		for(int i=0;i<n;i++) {for(int j=0;j<m;j++) {x.e[0][i]+=fabs(e[j][i]);}}
		for(int i=0;i<m;i++) {if(x.e[0][i]>s) s=x.e[0][i];}
		break;
	  case 'I' :
		for(int i=0;i<m;i++) {for(int j=0;j<n;j++) {x.e[0][i]+=fabs(e[i][j]);}}
		for(int i=0;i<m;i++) {if(x.e[0][i]>s) s=x.e[0][i];}
		break;
	  case '2' :
		//coming soon!
		break;
	  case 'F' :
		for(int i=0;i<m;i++) {for(int j=0;j<n;j++) {s+=pow(e[i][j],2);}}
		s=sqrt(s);
		break;
	  default :
		s=-1.;
		break;
    }
    return s;
}

matrix matrix::T()
{
    matrix t;
    t.init(n,m);
    for(int i=0;i<m;i++) {for(int j=0;j<n;j++) {t.e[j][i]=e[i][j];}}
    return t;
}

double matrix::Tr()
{
    if(m!=n) {cout<<"Matrix Trace Calculation Failure : Matrix Not Square."<<endl;exit(1);}
    double s=0.;
    for(int i=0;i<m;i++) {s+=e[i][i];}
    return s;
}

matrix matrix::reduce(int r, int c)
{
    matrix s;
    int a,b,h,k;
    a=m-1;
    b=n-1;
    if((a<=0)||(b<=0)) {s=(*this);cout<<"Matrix Cannot Be Reduced Further.";}
    else 
    {
	  s.init(a,b);
	  for(int i=0;i<m;i++) {for(int j=0;j<n;j++) if((i!=r)&&(j!=c))
	  {{
		h=i;
		k=j;
		if(i>r) {h=i-1;}
		if(j>c) {k=j-1;}
		s.e[h][k]=e[i][j];
	  }}}
    }
    return s;
}

double matrix::det()
{
    if(m!=n) {cout<<"Matrix Determinant Calculation Failure : Matrix Not Square."<<endl;exit(1);}
    double s=0.,dm=0.;
    if(m>2)
    {
	  for(int i=0;i<m;i++) {matrix mij;mij=reduce(0,i);dm=mij.det();s+=(pow(-1.,i)*e[0][i]*dm);}
    }
    else if(m==2)
    {
	  s=(e[0][0]*e[1][1])-(e[0][1]*e[1][0]);
    }
    else if(m==1)
    {
	  s=e[0][0];
    }
    else
    {
	  cout<<"Matrix Determinant Calculation Failure : Dim < 1."<<endl;
	  exit(1);
    }
    return s;
}

void matrix::Pivot(int a,int b)
{
    double t;
    for(int i=0;i<n;i++) {t=e[a][i];e[a][i]=e[b][i];e[b][i]=t;}
}

double matrix::fround(double f,int s)
{
    int n=1;
    double x,r,w1,w2;
    x=f;
    w1=pow(10.,s+1);
    w2=pow(10.,s-1);
    if(x!=0.)
    {
	  while (!((fabs(x)<w1)&&(fabs(x)>w2)))
	  {
		if(fabs(x)>w1) {r=0.1;} else {r=10.;}
		x*=r;
		n++;
	  }
	  n--;
	  x=floor(x+0.5);
	  r=pow(10.,n);
	  if(fabs(f)<w1) {x/=r;} else {x*=r;}
    }
    return x;
}

void matrix::FactorLU(matrix &U, matrix &L)
{
    if(det()==0.) {cout<<"Zero Determinant. Terminating Program."<<endl;exit(1);}
    L.init(m,n);
    U=(*this);
    for(int j=0;j<n;j++) {for (int i=j;i<m;i++)
    {
	  L.e[i][j]=U.e[i][j]/U.e[j][j];
	  if(i>j)
	  {
		for(int k=j;k<n;k++) {U.e[i][k]-=L.e[i][j]*U.e[j][k];}
	  }
    }}
}

void matrix::FactorLU(matrix &U, matrix &L,int s)
{
    if(det()==0.) {cout<<"Zero Determinant. Terminating Program."<<endl;exit(1);}
    L.init(m,n);
    U=(*this);
    double t;
    for(int j=0;j<n;j++)
    {
	  for (int i=j;i<m;i++)
	  {
		L.e[i][j]=U.e[i][j]/U.e[j][j];
		L.e[i][j]=fround(L.e[i][j],s);
		if(i>j)
		{
		    U.e[i][j]=0.;
		    for(int k=j+1;k<n;k++) {t=L.e[i][j]*U.e[j][k];t=fround(t,s);U.e[i][k]-=t;U.e[i][k]=fround(U.e[i][k],s);}
		}
	  }
    }
}

void matrix::FactorLU(matrix &U, matrix &L,matrix &P)
{
    if(det()==0.) {cout<<"Zero Determinant. Terminating Program."<<endl;exit(1);}
    L.init(m,n);
    P.init(m,n);
    P.ident();
    U=(*this);
    for(int j=0;j<n;j++)
    {
	  int r=j;
	  for(int i=j+1;i<m;i++) {if(fabs(U.e[i][j])>fabs(U.e[r][j])) {r=i;}}
	  if(r!=j) {L.Pivot(j,r);U.Pivot(j,r);P.Pivot(j,r);}
	  for (int i=j;i<m;i++)
	  {
		L.e[i][j]=U.e[i][j]/U.e[j][j];
		if(i>j)
		{
		    for(int k=j;k<n;k++) {U.e[i][k]-=L.e[i][j]*U.e[j][k];}
		}
	  }
    }
}

void matrix::FactorLU(matrix &U, matrix &L,matrix &P,int s)
{
    if(det()==0.) {cout<<"Zero Determinant. Terminating Program."<<endl;exit(1);}
    L.init(m,n);
    P.init(m,n);
    P.ident();
    U=(*this);
    double t;
    for(int j=0;j<n;j++)
    {
	  int r=j;
	  for(int i=j+1;i<m;i++) {if(fabs(U.e[i][j])>fabs(U.e[r][j])) {r=i;}}
	  if(r!=j) {L.Pivot(j,r);U.Pivot(j,r);P.Pivot(j,r);}
	  for (int i=j;i<m;i++)
	  {
		L.e[i][j]=U.e[i][j]/U.e[j][j];
		L.e[i][j]=fround(L.e[i][j],s);
		if(i>j)
		{
		    U.e[i][j]=0.;
		    for(int k=j+1;k<n;k++) {t=L.e[i][j]*U.e[j][k];t=fround(t,s);U.e[i][k]-=t;U.e[i][k]=fround(U.e[i][k],s);}
		}
	  }
    }
}

matrix matrix::GaussElimLU(matrix &U, matrix &L, matrix &b)
{
    if(((U.m==0)&&(U.n==0))&&((L.m==0)&&(L.n==0))) {FactorLU(U,L);}
    matrix x;
    x=b;
    for(int i=0;i<m;i++) {for(int j=0;j<i;j++) {x.e[i][0]-=L.e[i][j]*x.e[j][0];}}
    for(int i=m-1;i>=0;i--) {for(int j=i+1;j<n;j++) {x.e[i][0]-=U.e[i][j]*x.e[j][0];}x.e[i][0]/=U.e[i][i];}
    return x;
}

matrix matrix::GaussElimLU(matrix &U, matrix &L, matrix &b,int s)
{
    if(((U.m==0)&&(U.n==0))&&((L.m==0)&&(L.n==0))) {FactorLU(U,L,s);}
    matrix x;
    x=b;
    for(int i=0;i<m;i++) {for(int j=0;j<i;j++) {x.e[i][0]-=L.e[i][j]*x.e[j][0];x.e[i][0]=fround(x.e[i][0],s);}}
    for(int i=m-1;i>=0;i--) {for(int j=i+1;j<n;j++) {x.e[i][0]-=U.e[i][j]*x.e[j][0];x.e[i][0]=fround(x.e[i][0],s);}x.e[i][0]/=U.e[i][i];x.e[i][0]=fround(x.e[i][0],s);}
    return x;
}

matrix matrix::GaussElimLU(matrix &U, matrix &L, matrix &P, matrix &b)
{
    if(((U.m==0)&&(U.n==0))&&((L.m==0)&&(L.n==0))) {FactorLU(U,L,P);}
    matrix x;
    x=P*b;
    for(int i=0;i<m;i++) {for(int j=0;j<i;j++) {x.e[i][0]-=L.e[i][j]*x.e[j][0];}}
    for(int i=m-1;i>=0;i--) {for(int j=i+1;j<n;j++) {x.e[i][0]-=U.e[i][j]*x.e[j][0];}x.e[i][0]/=U.e[i][i];}
    return x;
}

matrix matrix::GaussElimLU(matrix &U, matrix &L, matrix &P, matrix &b,int s)
{
    if(((U.m==0)&&(U.n==0))&&((L.m==0)&&(L.n==0))) {FactorLU(U,L,P,s);}
    matrix x;
    x=P*b;
    for(int i=0;i<m;i++) {for(int j=0;j<i;j++) {x.e[i][0]-=L.e[i][j]*x.e[j][0];x.e[i][0]=fround(x.e[i][0],s);}}
    for(int i=m-1;i>=0;i--) {for(int j=i+1;j<n;j++) {x.e[i][0]-=U.e[i][j]*x.e[j][0];x.e[i][0]=fround(x.e[i][0],s);}x.e[i][0]/=U.e[i][i];x.e[i][0]=fround(x.e[i][0],s);}
    return x;
}

bool matrix::Symmetric() {return(T()==*this);}

bool matrix::PosDef()
{
    if(m!=n) {cout<<"Positive Definite Check Failure : Matrix Not Square."<<endl;exit(1);}
    bool c=true;
    for(int i=1;i<=m;i++)
    {
	  matrix ai;
	  ai.init(i,i);
	  for(int j=0;j<i;j++) {for(int k=0;k<i;k++) {ai.e[j][k]=e[j][k];}}
	  if(ai.det()<=0) {c=false;break;}
    }
    return c;
}

matrix matrix::CholeskyL()
{
    if(!Symmetric()) {cout<<"Cholesky Calculation Failure : Matrix Not Symmetric."<<endl;exit(1);}
    if(!PosDef()) {cout<<"Cholesky Calculation Failure : Matrix Not Positive Definite."<<endl;exit(1);}
    matrix L;
    L.init(m,m);
    for(int i=0;i<m;i++)
    {
	  for(int j=0;j<i;j++)
	  {
		for(int k=0;k<j;k++) {L.e[i][j]+=L.e[i][k]*L.e[j][k];}L.e[i][j]=e[i][j]-L.e[i][j];L.e[i][j]/=L.e[j][j];
	  }
	  for(int k=0;k<i;k++) {L.e[i][i]+=pow(L.e[i][k],2);}L.e[i][i]=e[i][i]-L.e[i][i];L.e[i][i]=sqrt(L.e[i][i]);
    }
    return L;
}

bool matrix::Tridiagonal()
{
    if(m!=n) {cout<<"Tridiagonal Check Failure : Matrix Not Square."<<endl;exit(1);}
    bool c=true;
    if(m>2)
    {
	  for(int i=2;i<m;i++) {for(int j=0;j<i-1;j++) {if((e[i][j]!=0.)||(e[j][i]!=0.)) {c=false;}}}
    }
    return c;
}

void matrix::FactorTLU(matrix &U,matrix &L)
{
    if(det()==0.) {cout<<"Zero Determinant. Terminating Program."<<endl;exit(1);}
    if(!Tridiagonal()) {cout<<"Tridiagonal Factorization Failure : Matrix Not Tridiagonal."<<endl;exit(1);}
    L.init(m,n);
    U.init(m,n);
    U.ident();
    for(int i=0;i<m;i++) {L.e[i][i]=e[i][i];if(i>0) {L.e[i][i-1]=e[i][i-1];L.e[i][i]-=(L.e[i][i-1]*U.e[i-1][i]);}if(i<m-1) {U.e[i][i+1]=e[i][i+1]/L.e[i][i];}}
}

matrix matrix::TridiagonalSolve(matrix &U, matrix &L, matrix &b)
{
    if(((U.m==0)&&(U.n==0))&&((L.m==0)&&(L.n==0))) {FactorTLU(U,L);}
    matrix x;
    x=b;
    for(int i=0;i<m;i++) {if(i>0) {x.e[i][0]-=(L.e[i][i-1]*x.e[i-1][0]);}x.e[i][0]/=L.e[i][i];}
    for(int i=m-2;i>=0;i--) {x.e[i][0]-=(U.e[i][i+1]*x.e[i+1][0]);}
    return x;
}
