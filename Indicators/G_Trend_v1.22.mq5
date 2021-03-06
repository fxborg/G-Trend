//+------------------------------------------------------------------+
//|                                               G_Trend_v1.22.mq5  |
//| G_Trend_v1.22                             Copyright 2016, fxborg |
//|                                   http://fxborg-labo.hateblo.jp/ |
//|This indicator drawing trend lines automatically using centroid of|
//|convex hulls.                                                     |
//|http://fxborg-labo.hateblo.jp/archive/category/Auto%20Trend%20Line|
//+------------------------------------------------------------------+
#property copyright "Copyright 2016, fxborg"
#property link      "http://fxborg-labo.hateblo.jp/"
#property version   "1.22"
#property indicator_chart_window


#property indicator_buffers 7
#property indicator_plots 1

#property indicator_type1         DRAW_LINE 
#property indicator_color1        clrRed
#property indicator_width1 2

input int InpOID=1; //  Object ID
input int InpConvexPeriod=40; //  Polygon Period
input int InpRegrPeriod=8;    //  Regression Period
input int Inp1stPeriod=60;    //  Trend Period
input int  InpShowLine=1;    //Show Line (1:show ,0:hide)  
input color Inp1stColor=clrRed; // Color 1
input color Inp2ndColor=clrMediumSlateBlue; // Color 2
input color Inp3rdColor=clrAqua; // Color 3
input int InpMaxBars=50000; // MaxBars

double HI[];
double LO[];
double CX[];
double CY[];
double LA[];
double LB[];
double TREND[];

int WinNo=ChartWindowFind();
int min_rates_total=InpConvexPeriod+Inp1stPeriod;
//+------------------------------------------------------------------+
//| Custom indicator initialization function                         |
//+------------------------------------------------------------------+
int OnInit()
  {
   ObjectDeleteByName(StringFormat("GTrend_%d",InpOID));
   if(InpShowLine==0)
      PlotIndexSetInteger(0,PLOT_DRAW_TYPE,DRAW_NONE);
   else
      PlotIndexSetInteger(0,PLOT_DRAW_TYPE,DRAW_LINE);

   SetIndexBuffer(0,TREND,INDICATOR_DATA);//---
   SetIndexBuffer(1,CX,INDICATOR_CALCULATIONS);
   SetIndexBuffer(2,CY,INDICATOR_CALCULATIONS);
   SetIndexBuffer(3,LA,INDICATOR_CALCULATIONS);
   SetIndexBuffer(4,LB,INDICATOR_CALCULATIONS);
   SetIndexBuffer(5,HI,INDICATOR_CALCULATIONS);
   SetIndexBuffer(6,LO,INDICATOR_CALCULATIONS);

///  --- 
//--- digits
   IndicatorSetInteger(INDICATOR_DIGITS,2);
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
   for(i=first; i<rates_total-1 && !IsStopped(); i++)
     {
      CX[i]=EMPTY_VALUE;
      CY[i]=EMPTY_VALUE;
      LA[i]=EMPTY_VALUE;
      LB[i]=EMPTY_VALUE;
      TREND[i]=EMPTY_VALUE;

      //---
      HI[i-1]=(high[i-2]+high[i-1]+high[i])/3;
      LO[i-1]=(low[i-2]+low[i-1]+low[i])/3;
      //---
      if(i<rates_total-InpMaxBars)continue;
      //---
      int i1st=begin_pos+InpConvexPeriod*2;
      if(i<=i1st)continue;
      //---

      double upper[][2];
      double lower[][2];

      //---
      convex_hull(upper,lower,HI,LO,i-1,InpConvexPeriod);
      //---

      int up_sz=int(ArraySize(upper)*0.5);
      int lo_sz=int(ArraySize(lower)*0.5);
      calc_vector(upper,lower,i-1);

      //---
      int i2nd=i1st+Inp1stPeriod*2;
      if(i<=i2nd)continue;
      //---
      double alpha,y0,y1,y2;
      int from_x,x1;
      //---
      calc_trend(alpha,from_x,y0,x1,y1,y2,time,i-1);
      //---
      //--- Trend
      int x1_len=i-x1;
      for(int j=0;j<=x1_len+1;j++) TREND[x1+j]=y1+(alpha*j);
      drawTrend(InpOID,1,Inp1stColor,from_x,y0,i-1,y2,time,STYLE_DOT,1);

      //---
      double upper2[][2];
      double lower2[][2];
      int up_n=0;
      int lo_n=0;
      copy_outside(upper2,up_n,lower2,lo_n,HI,LO,alpha,from_x,y0,i);

      //---
      double ra=0,rb=0;
      double r_ma=0;
      double r_max=0;
      calc_snr(ra,rb,r_ma,r_max,upper2,up_n,i,1);
      //---
      double sa=0,sb=0;
      double s_ma=0;
      double s_max=0;
      calc_snr(sa,sb,s_ma,s_max,lower2,lo_n,i,-1);

      double ry0,sy0,ry2,sy2;

      ry0=ra*from_x+rb+r_ma;
      ry2=ra*(i-1)+rb+r_ma;

      sy0=sa*from_x+sb-s_ma;
      sy2=sa*(i-1)+sb-s_ma;
      //---
      drawTrend(InpOID,2,Inp2ndColor,from_x,ry0,i-1,ry2,time,STYLE_DOT,1);
      drawTrend(InpOID,3,Inp2ndColor,from_x,sy0,i-1,sy2,time,STYLE_DOT,1);

      //---
      ry0=ra*from_x+rb+r_max;
      ry2=ra*(i-1)+rb+r_max;

      sy0=sa*from_x+sb-s_max;
      sy2=sa*(i-1)+sb-s_max;
      //---
      drawTrend(InpOID,4,Inp3rdColor,from_x,ry0,i-1,ry2,time,STYLE_SOLID,2);
      drawTrend(InpOID,5,Inp3rdColor,from_x,sy0,i-1,sy2,time,STYLE_SOLID,2);
      //---

     }

//--- return value of prev_calculated for next call
   return(rates_total);
  }
//---
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
void calc_snr(double &a,double &b,double &ma,double &max,const double &points[][2],const int sz,const int i,const int dir)
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
void calc_trend(double &alpha,int &x0,double &y0,int &x1,double &y1,double &y2,const datetime  &time[],const int i)
  {
   double sumx=0;
   double sumy=0;
   double a=0;
   int a_count=0;
   int ifrom=0;
   int cnt=0;
   int len=Inp1stPeriod;
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
   int from_x=int(CX[ifrom]-InpConvexPeriod*0.5);
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
void calc_vector(double  &upper[][2],double  &lower[][2],const int i)
  {
//if(CX[i]!=EMPTY_VALUE)return;
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
void drawTrend(const int oid,int no,
               const color clr,const int x0,const double y0,const int x1,const double y1,
               const datetime &time[],const ENUM_LINE_STYLE style,const int width)
  {

   if(-1<ObjectFind(0,StringFormat("GTrend_%d_%d",oid,no)))
     {
      ObjectMove(0,StringFormat("GTrend_%d_%d",oid,no),0,time[x0],y0);
      ObjectMove(0,StringFormat("GTrend_%d_%d",oid,no),1,time[x1],y1);
     }
   else
     {
      ObjectCreate(0,StringFormat("GTrend_%d_%d",oid,no),OBJ_TREND,WinNo,time[x0],y0,time[x1],y1);
      ObjectSetInteger(0,StringFormat("GTrend_%d_%d",oid,no),OBJPROP_COLOR,clr);
      ObjectSetInteger(0,StringFormat("GTrend_%d_%d",oid,no),OBJPROP_STYLE,style);
      ObjectSetInteger(0,StringFormat("GTrend_%d_%d",oid,no),OBJPROP_WIDTH,width);
      ObjectSetInteger(0,StringFormat("GTrend_%d_%d",oid,no),OBJPROP_RAY_RIGHT,true);
     }
  }
//+------------------------------------------------------------------+
//|                                                                  |
//+------------------------------------------------------------------+
void ObjectDeleteByName(string prefix)
  {
   int total=ObjectsTotal(0),
   length=StringLen(prefix);
   for(int i=total-1; i>=0; i--)
     {
      string objName=ObjectName(0,i);
      if(StringSubstr(objName,0,length)==prefix)
        {
         ObjectDelete(0,objName);
        }
     }
  }
//+------------------------------------------------------------------+
