//+------------------------------------------------------------------+
//|                                           G_Trend_deb_v1.00.mq5  |
//| G_Trend_deb_v1.00                         Copyright 2016, fxborg |
//|                                   http://fxborg-labo.hateblo.jp/ |
//|This indicator drawing trend lines automatically using centroid of|
//|convex hulls.                                                     |
//|http://fxborg-labo.hateblo.jp/archive/category/Auto%20Trend%20Line|
//+------------------------------------------------------------------+
#property copyright "Copyright 2016, fxborg"
#property link      "http://fxborg-labo.hateblo.jp/"
#property version   "1.00"
#property indicator_chart_window


#property indicator_buffers 23
#property indicator_plots 10

#property indicator_type1         DRAW_LINE
#property indicator_color1       clrRed
#property indicator_width1 1


#property indicator_type2     DRAW_LINE 
#property indicator_color2    clrGold
#property indicator_width2 2

#property indicator_type3         DRAW_LINE 
#property indicator_color3       clrSilver
#property indicator_width3    1
#property indicator_style3    STYLE_DOT

#property indicator_type4         DRAW_LINE 
#property indicator_color4       clrGold
#property indicator_width4 2

#property indicator_type5         DRAW_LINE 
#property indicator_color5       clrSilver
#property indicator_width5    1
#property indicator_style5    STYLE_DOT

#property indicator_type6         DRAW_LINE 
#property indicator_color6       clrCyan
#property indicator_width6 2

#property indicator_type7         DRAW_LINE 
#property indicator_color7       clrSilver
#property indicator_width7    1
#property indicator_style7    STYLE_DOT

#property indicator_type8         DRAW_LINE 
#property indicator_color8       clrCyan
#property indicator_width8 2

#property indicator_type9         DRAW_LINE 
#property indicator_color9       clrSilver
#property indicator_width9    1
#property indicator_style9    STYLE_DOT

#property indicator_type10         DRAW_LINE
#property indicator_color10       clrMagenta
#property indicator_width10 1



input int Inp1stConvexPeriod=40; //  Polygon Period
input int Inp1stTrendPeriod=50;    //  Trend Period
input int Inp2ndConvexPeriod=60; //  Polygon Period
input int Inp2ndTrendPeriod=200;    //  Trend Period

input int InpMaxBars=20000; // MaxBars
input int InpRegrPeriod=8;    //  Regression Period

double SIG[];
double SIG_CLR[];

double HI[];
double LO[];
double CX1[];
double CY1[];
double LA1[];
double LB1[];
double TREND1[];
double CX2[];
double CY2[];
double LA2[];
double LB2[];
double TREND2[];

double FR1[];
double FR2[];
double FS1[];
double FS2[];
double SR1[];
double SR2[];
double SS1[];
double SS2[];
int min_rates_total=MathMax(Inp1stConvexPeriod+Inp1stTrendPeriod,Inp2ndConvexPeriod+Inp2ndTrendPeriod);
int Span1=Inp1stTrendPeriod+int(Inp1stConvexPeriod/2);
int Span2=Inp2ndTrendPeriod+int(Inp2ndConvexPeriod/2);
//+------------------------------------------------------------------+
//| Custom indicator initialization function                         |
//+------------------------------------------------------------------+
int OnInit()
  {

   int i=0;
   SetIndexBuffer(i++,TREND1,INDICATOR_DATA);//---
   SetIndexBuffer(i++,FR2,INDICATOR_DATA);//---
   SetIndexBuffer(i++,FR1,INDICATOR_DATA);//---
   SetIndexBuffer(i++,FS2,INDICATOR_DATA);//---
   SetIndexBuffer(i++,FS1,INDICATOR_DATA);//---
   SetIndexBuffer(i++,SR2,INDICATOR_DATA);//---
   SetIndexBuffer(i++,SR1,INDICATOR_DATA);//---
   SetIndexBuffer(i++,SS2,INDICATOR_DATA);//---
   SetIndexBuffer(i++,SS1,INDICATOR_DATA);//---
   SetIndexBuffer(i++,TREND2,INDICATOR_DATA);//---

   SetIndexBuffer(i++,CX1,INDICATOR_CALCULATIONS);
   SetIndexBuffer(i++,CX2,INDICATOR_CALCULATIONS);
   SetIndexBuffer(i++,CY1,INDICATOR_CALCULATIONS);
   SetIndexBuffer(i++,CY2,INDICATOR_CALCULATIONS);
   SetIndexBuffer(i++,LA1,INDICATOR_CALCULATIONS);
   SetIndexBuffer(i++,LA2,INDICATOR_CALCULATIONS);
   SetIndexBuffer(i++,LB1,INDICATOR_CALCULATIONS);
   SetIndexBuffer(i++,LB2,INDICATOR_CALCULATIONS);
   SetIndexBuffer(i++,HI,INDICATOR_CALCULATIONS);
   SetIndexBuffer(i++,LO,INDICATOR_CALCULATIONS);

///  --- 
//--- digits
   IndicatorSetInteger(INDICATOR_DIGITS,3);
   return(0);
  }
//+------------------------------------------------------------------+
//| Custom indicator iteration function                              |
//+------------------------------------------------------------------+
int OnCalculate(const int rates_total,
                const int prev_calculated,
                const datetime &time[],
                const double &open[],
                const double &high[],
                const double &low[],
                const double &close[],
                const long &tick_volume[],
                const long &volume[],
                const int &spread[])
  {
//---
   int i,first;
   if(rates_total<=min_rates_total) return(0);
//---
   int begin_pos=min_rates_total;

//---
   first=begin_pos;
   if(first+1<prev_calculated) first=prev_calculated-2;


//---
   for(i=first; i<rates_total && !IsStopped(); i++)
     {
      SR1[i-Span2-1]=EMPTY_VALUE;
      SR2[i-Span2-1]=EMPTY_VALUE;
      SS1[i-Span2-1]=EMPTY_VALUE;
      SS2[i-Span2-1]=EMPTY_VALUE;
      FR1[i-Span1-1]=EMPTY_VALUE;
      FR2[i-Span1-1]=EMPTY_VALUE;
      FS1[i-Span1-1]=EMPTY_VALUE;
      FS2[i-Span1-1]=EMPTY_VALUE;
      SR1[i]=EMPTY_VALUE;
      SR2[i]=EMPTY_VALUE;
      SS1[i]=EMPTY_VALUE;
      SS2[i]=EMPTY_VALUE;
      FR1[i]=EMPTY_VALUE;
      FR2[i]=EMPTY_VALUE;
      FS1[i]=EMPTY_VALUE;
      FS2[i]=EMPTY_VALUE;


      CX1[i]=EMPTY_VALUE;
      CY1[i]=EMPTY_VALUE;
      LA1[i]=EMPTY_VALUE;
      LB1[i]=EMPTY_VALUE;
      TREND1[i]=EMPTY_VALUE;
      CX2[i]=EMPTY_VALUE;
      CY2[i]=EMPTY_VALUE;
      LA2[i]=EMPTY_VALUE;
      LB2[i]=EMPTY_VALUE;
      TREND2[i]=EMPTY_VALUE;

      //---
      HI[i-1]=(high[i-2]+high[i-1]+high[i])/3;
      LO[i-1]=(low[i-2]+low[i-1]+low[i])/3;
      //---
      if(i<rates_total-InpMaxBars)continue;
      //---
      int i1st=begin_pos+MathMax(Inp1stConvexPeriod,Inp2ndConvexPeriod)*2;
      if(i<=i1st)continue;
      //---

      double upper1[][2];
      double lower1[][2];
      double upper2[][2];
      double lower2[][2];

      //---
      convex_hull(upper1,lower1,HI,LO,i-1,Inp1stConvexPeriod);
      convex_hull(upper2,lower2,HI,LO,i-1,Inp2ndConvexPeriod);
      //---

      int up1_sz=int(ArraySize(upper1)*0.5);
      int lo1_sz=int(ArraySize(lower1)*0.5);
      calc_vector(CX1,CY1,LA1,LB1,upper1,lower1,i-1);
      //---
      int up2_sz=int(ArraySize(upper2)*0.5);
      int lo2_sz=int(ArraySize(lower2)*0.5);
      calc_vector(CX2,CY2,LA2,LB2,upper2,lower2,i-1);
      //---
      int i2nd=i1st+MathMax(Inp1stTrendPeriod,Inp2ndTrendPeriod)*2;
      if(i<=i2nd)continue;
      //---

      //---
      double alpha_1,y0_1,y1_1,y2_1;
      int from_x_1,x1_1;
      double alpha_2,y0_2,y1_2,y2_2;
      int from_x_2,x1_2;
      //---
      calc_trend(alpha_1,from_x_1,y0_1,x1_1,y1_1,y2_1,CX1,CY1,LA1,time,i-1,Inp1stConvexPeriod,Inp1stTrendPeriod);
      calc_trend(alpha_2,from_x_2,y0_2,x1_2,y1_2,y2_2,CX2,CY2,LA2,time,i-1,Inp2ndConvexPeriod,Inp2ndTrendPeriod);
      //---
      //--- Trend
      for(int j=0;j<=(i-x1_1);j++) TREND1[x1_1+j]=y1_1+(alpha_1*j);
      for(int j=0;j<=(i-x1_2);j++) TREND2[x1_2+j]=y1_2+(alpha_2*j);


      //---
      double ra_1=0,rb_1=0,r1_1=0,r2_1=0;
      double sa_1=0,sb_1=0,s1_1=0,s2_1=0;

      //---
      calc_snr(ra_1,rb_1,r1_1,r2_1,sa_1,sb_1,s1_1,s2_1,HI,LO,alpha_1,from_x_1,y0_1,i);
      //---

      //---

      double ra_2=0,rb_2=0,r1_2=0,r2_2=0;
      double sa_2=0,sb_2=0,s1_2=0,s2_2=0;

      //---
      calc_snr(ra_2,rb_2,r1_2,r2_2,sa_2,sb_2,s1_2,s2_2,HI,LO,alpha_2,from_x_2,y0_2,i);
      //---

      double _r1_0= ra_1 * (from_x_1) + rb_1 + r1_1;
      double _s1_0= sa_1 * (from_x_1) + sb_1 - s1_1;

      double _r1_1= ra_1 * (i-1) + rb_1 + r1_1;
      double _r2_1= ra_1 * (i-1) + rb_1 + r2_1;
      double _s1_1= sa_1 * (i-1) + sb_1 - s1_1;
      double _s2_1= sa_1 * (i-1) + sb_1 - s2_1;

      double _r1_2= ra_2 * (i-1) + rb_2 + r1_2;
      double _r2_2= ra_2 * (i-1) + rb_2 + r2_2;
      double _s1_2= sa_2 * (i-1) + sb_2 - s1_2;
      double _s2_2= sa_2 * (i-1) + sb_2 - s2_2;

      //--- Trend
      for(int j=i-Span1;j<i;j++) FR1[j]=(ra_1*j)+rb_1+r1_1;
      for(int j=i-Span1;j<i;j++) FR2[j]=(ra_1*j)+rb_1+r2_1;
      for(int j=i-Span1;j<i;j++) FS1[j]=(sa_1*j)+sb_1-s1_1;
      for(int j=i-Span1;j<i;j++) FS2[j]=(sa_1*j)+sb_1-s2_1;

      for(int j=i-Span2;j<i;j++) SR1[j]=(ra_2*j)+rb_2+r1_2;
      for(int j=i-Span2;j<i;j++) SR2[j]=(ra_2*j)+rb_2+r2_2;
      for(int j=i-Span2;j<i;j++) SS1[j]=(sa_2*j)+sb_2-s1_2;
      for(int j=i-Span2;j<i;j++) SS2[j]=(sa_2*j)+sb_2-s2_2;


     }

//--- return value of prev_calculated for next call
   return(rates_total);
  }
//---
void calc_snr(double &ra,double &rb,double &r1,double &r2,
              double &sa,double &sb,double &s1,double &s2,
              const double  &hi[],const double  &lo[],
              const double alpha,const int from_x,const double y0,const int i)
  {
//---
   double upper[][2];
   double lower[][2];
   int up_n=0;
   int lo_n=0;
   copy_outside(upper,up_n,lower,lo_n,hi,lo,alpha,from_x,y0,i);

//---
   _calc_snr(ra,rb,r1,r2,upper,up_n,i,1);
//---
   _calc_snr(sa,sb,s1,s2,lower,lo_n,i,-1);

  }
//+------------------------------------------------------------------+
//|                                                                  |
//+------------------------------------------------------------------+
void copy_outside(double &upper[][2],int &up_n,double &lower[][2],int &lo_n,const double &hi[],const double &lo[],const double alpha,const int from_x,const double y0,const int i)
  {
   int sz=i-from_x;
   ArrayResize(upper,sz);
   ArrayResize(lower,sz);

   double trend=y0;
   up_n=0;
   lo_n=0;
   for(int j=from_x;j<=i-1;j++)
     {

      if(hi[j]>trend)
        {
         upper[up_n][0]=j;
         upper[up_n][1]=hi[j];
         up_n++;
        }
      if(lo[j]<trend)
        {
         lower[lo_n][0]=j;
         lower[lo_n][1]=lo[j];
         lo_n++;
        }
      trend+=alpha;
     }
  }
//+------------------------------------------------------------------+
//|                                                                  |
//+------------------------------------------------------------------+
void _calc_snr(double &a,double &b,double &ma,double &max,const double &points[][2],const int sz,const int i,const int dir)
  {
//---
   _regression(a,b,points,sz);
   ma=0;
   max=0;
   int cnt=0;
//---
   for(int j=0;j<sz;j++)
     {
      double v=dir *(points[j][1]-(a*points[j][0]+b));
      if(v>=0)
        {
         if(v>max)max=v;
         ma+=v;
         cnt++;
        }
     }
   ma/=MathMax(1,cnt);

   if(ma*2.5<max)max=ma*2.5;

  }
//+------------------------------------------------------------------+
//|                                                                  |
//+------------------------------------------------------------------+
bool  get_slope(double &slope,double &y0,double &y1,const double &vertexs[][2],const double alpha,const int from_x,const int i)
  {
   int sz=int(ArraySize(vertexs)*0.5);
//---
   double df=999999999999;
   int pos=-1;
//---
   for(int j=0;j<sz-1;j++)
     {
      double span=vertexs[j][0]-vertexs[j+1][0];
      if(span<=0)continue;
      double tmp=(vertexs[j][1]-vertexs[j+1][1])/span;
      if(df>MathAbs(alpha-tmp))
        {
         df=MathAbs(alpha-tmp);
         pos=j;
         slope=tmp;
        }
     }
//---
   if(pos == -1)return false;
   y1=(i-vertexs[pos][0])*slope+vertexs[pos][1];
   y0=y1-(i-from_x)*slope;
   return true;
  }
//+------------------------------------------------------------------+
//|                                                                  |
//+------------------------------------------------------------------+
void convex_hull(double &upper[][2],double &lower[][2],const  double &high[],const double &low[],const int i,const int len)
  {

   ArrayResize(upper,len,len);
   int k=0;
   for(int j=0;j<len;j++)
     {
      while(k>=2 && 
            (cross(upper[k-2][0],upper[k-2][1],
            upper[k-1][0],upper[k-1][1],
            i-j,high[i-j]))<=0)
        {
         k--;
        }

      upper[k][0]= i-j;
      upper[k][1]= high[i-j];
      k++;
     }
   ArrayResize(upper,k,len);

   ArrayResize(lower,len,len);
   k=0;
   for(int j=0;j<len;j++)
     {
      while(k>=2 && 
            (cross(lower[k-2][0],lower[k-2][1],
            lower[k-1][0],lower[k-1][1],
            i-j,low[i-j]))>=0)
        {
         k--;
        }

      lower[k][0]= i-j;
      lower[k][1]= low[i-j];
      k++;
     }
   ArrayResize(lower,k,len);

  }
//+------------------------------------------------------------------+
//|                                                                  |
//+------------------------------------------------------------------+
void calc_trend(double &alpha,int &x0,double &y0,int &x1,double &y1,double &y2
                ,const double  &CX[],const double  &CY[],const double  &LA[]
                ,const datetime  &time[],const int i
                ,const int convex_period,const int period)
  {
   double sumx=0;
   double sumy=0;
   double a=0;
   int a_count=0;
   int ifrom=0;
   int cnt=0;
   int len=period;
   for(int j=0;j<=len;j++)
     {

      if(CX[i-j]!=EMPTY_VALUE && LA[i-j]!=EMPTY_VALUE)
        {
         a+=LA[i-j];
         ifrom=i-j;
         a_count++;
         sumx+=CX[ifrom];
         sumy+=CY[ifrom];

        }
     }
   double ax=i-(sumx/a_count);
   double ay=sumy/a_count;
   int from_x=int(CX[ifrom]-convex_period*0.5);
   double aa=(a/a_count);
   double y=aa*ax+ay;
   double span=i-from_x;
   double from_y=y-aa*span;
   alpha=aa;
   x0=from_x;
   y0=from_y;
   x1=int((sumx/a_count)+0.5);
   y1=y-aa*(i-x1);
   y2=y;



  }
//+------------------------------------------------------------------+
//|                                                                  |
//+------------------------------------------------------------------+
void calc_vector(double  &CX[],double  &CY[],double  &LA[],double  &LB[],
                 double  &upper[][2],double  &lower[][2],const int i)
  {
   if(CX[i]!=EMPTY_VALUE)return;
   int up_sz=int(ArraySize(upper)*0.5);
   int lo_sz=int(ArraySize(lower)*0.5);


   double mx,my;
   calc_centroid(mx,my,upper,lower);
   if(mx<i)
     {

      CY[i]=my;
      CX[i]=mx;
      double a,b;
      regression(a,b,CX,CY,i-InpRegrPeriod-1,i);
      LA[i]=a;
      LB[i]=b;
     }
  }
//+------------------------------------------------------------------+
//|                                                                  |
//+------------------------------------------------------------------+
void calc_centroid(double  &x,double  &y,const double  &upper[][2],const double  &lower[][2])
  {
   double vertices[][2];
   int up_sz=int(ArraySize(upper)*0.5);
   int lo_sz=int(ArraySize(lower)*0.5);
   int sz=up_sz+lo_sz;
   ArrayResize(vertices,sz,sz);
   int n=0;
   for(int j=up_sz-1;j>=0;j--)
     {
      vertices[n][0]=upper[j][0];
      vertices[n][1]=upper[j][1];
      n++;
     }

   for(int j=0;j<lo_sz;j++)
     {
      vertices[n][0]=lower[j][0];
      vertices[n][1]=lower[j][1];
      n++;
     }
   ArrayResize(vertices,n,sz);

   int v_cnt=n;
   y=0;
   x=0;
   double signedArea=0.0;
   double x0 = 0.0; // Current vertex X
   double y0 = 0.0; // Current vertex Y
   double x1 = 0.0; // Next vertex X
   double y1 = 0.0; // Next vertex Y
   double a = 0.0;  // Partial signed area

                    // For all vertices
   int i=0;
   for(i=0; i<v_cnt-1; i++)
     {
      x0 = vertices[i][0];
      y0 = vertices[i][1];
      if(i==v_cnt-2)
        {
         x1 = vertices[0][0];
         y1 = vertices[0][1];
        }
      else
        {
         x1 = vertices[i+1][0];
         y1 = vertices[i+1][1];
        }
      a=x0*y1-x1*y0;
      signedArea+=a;
      x += (x0 + x1)*a;
      y += (y0 + y1)*a;
     }
   if(signedArea!=0.0)
     {
      signedArea*=0.5;
      x /= (6.0*signedArea);
      y /= (6.0*signedArea);
     }
  }
//+------------------------------------------------------------------+
//|                                                                  |
//+------------------------------------------------------------------+
void regression(double  &a,double  &b,const double &x[],const double &y[],const int from,const int to)
  {
   int temp_sz=to-from;
   double temp[][2];
   ArrayResize(temp,temp_sz+1);
   int n=0;
   for(int k=from;k<=to;k++)
     {
      if(x[k]==EMPTY_VALUE)continue;
      if(y[k]==EMPTY_VALUE)continue;
      temp[n][0]=x[k];
      temp[n][1]=y[k];
      n++;
     }
   _regression(a,b,temp,n);
  }
//+------------------------------------------------------------------+
//|
//+------------------------------------------------------------------+
void _regression(double  &a,double  &b,const double &data[][2],const int cnt)
  {

   if(cnt==0)
     {
      a=EMPTY_VALUE;
      b=EMPTY_VALUE;
      return;
     }
//--- 
   double sumy=0.0; double sumx=0.0;
   double sumxy=0.0; double sumx2=0.0;

//--- 
   for(int n=0; n<cnt; n++)
     {
      //---
      sumx+=data[n][0];
      sumx2+= data[n][0]*data[n][0];
      sumy += data[n][1];
      sumxy+= data[n][0]*data[n][1];

     }
//---
   double c=sumx2-sumx*sumx/cnt;
   if(c==0.0)
     {
      a=0.0;
      b=sumy/cnt;
     }
   else
     {
      a=(sumxy-sumx*sumy/cnt)/c;
      b=(sumy-sumx*a)/cnt;
     }
  }
//+------------------------------------------------------------------+
//|                                                                  |
//+------------------------------------------------------------------+
double cross(const double ox,double oy,
             const double ax,double ay,
             const double bx,double by)
  {
   return ((ax - ox) * (by - oy) - (ay - oy) * (bx - ox));
  }
//+------------------------------------------------------------------+
//|                                                                  |
//+------------------------------------------------------------------+
