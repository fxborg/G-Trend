//+------------------------------------------------------------------+
//|                                                  AdaptiveTL.mq5  |
//| Adaptive TL                               Copyright 2016, fxborg |
//|                                   http://fxborg-labo.hateblo.jp/ |
//+------------------------------------------------------------------+
#property copyright "Copyright 2016, fxborg"
#property link      "http://fxborg-labo.hateblo.jp/"
#property version   "1.00"
#property indicator_chart_window


#property indicator_buffers 11
#property indicator_plots 4

#property indicator_type1         DRAW_LINE 
#property indicator_color1        clrRed
#property indicator_width1 2
#property indicator_type2         DRAW_LINE 
#property indicator_color2        clrDodgerBlue
#property indicator_width2 2

#property indicator_type3         DRAW_LINE 
#property indicator_color3        clrMediumSlateBlue
#property indicator_width3 1
#property indicator_style3       STYLE_DOT

#property indicator_type4         DRAW_LINE 
#property indicator_color4        clrMediumSlateBlue
#property indicator_width4 1
#property indicator_style4        STYLE_DOT

input int InpMaxPeriod=200;    //  MA  Period
input int InpPeriod=20;    //  Trend Period
double ATR[];
double HI[];
double LO[];
double S1[];
double R1[];
double S2[];
double R2[];
double VOLAT[];
int min_rates_total=InpPeriod+2;
//+------------------------------------------------------------------+
//| Custom indicator initialization function                         |
//+------------------------------------------------------------------+
int OnInit()
  {

   SetIndexBuffer(0,R2,INDICATOR_DATA);//---
   SetIndexBuffer(1,S2,INDICATOR_DATA);//---
   SetIndexBuffer(2,R1,INDICATOR_DATA);//---
   SetIndexBuffer(3,S1,INDICATOR_DATA);//---
   SetIndexBuffer(4,HI,INDICATOR_CALCULATIONS);
   SetIndexBuffer(5,LO,INDICATOR_CALCULATIONS);
   SetIndexBuffer(6,ATR,INDICATOR_CALCULATIONS);//---
   SetIndexBuffer(7,VOLAT,INDICATOR_DATA);

///  --- 
//--- digits
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
      ATR[i]=MathMax(high[i],close[i-1])-MathMin(low[i],close[i-1]);
      VOLAT[i]=fabs(close[i]-close[i-1]);
      S1[i]=EMPTY_VALUE;
      R1[i]=EMPTY_VALUE;
      S2[i]=EMPTY_VALUE;
      R2[i]=EMPTY_VALUE;
      if(i==rates_total-1)
        {
         S1[i]=S1[i-1]+(S1[i-1]-S1[i-2]);
         R1[i]=R1[i-1]+(R1[i-1]-R1[i-2]);
         S2[i]=S2[i-1]+(S2[i-1]-S2[i-2]);
         R2[i]=R2[i-1]+(R2[i-1]-R2[i-2]);
         continue;
        }
      if(i==begin_pos)continue;

      ATR[i]=(1-0.99) *ATR[i]+0.99*ATR[i-1];
      //---    
      if(i<=begin_pos+200)continue;

      //---
      HI[i-1]=(high[i-2]+high[i-1]+high[i])/3;
      LO[i-1]=(low[i-2]+low[i-1]+low[i])/3;


      //---
      int i1st=begin_pos+InpMaxPeriod+2;
      if(i<=i1st)continue;

      double limit=ATR[i]*InpPeriod;

      double dsum=0;
      double dmax=ATR[i]*4;
      double dmin=ATR[i]*0.5;

      int bk=1;

      for(bk=1;bk<InpMaxPeriod;bk++)
        {
         if(VOLAT[i-bk]<dmin)continue;
         dsum+=fmin(dmax,VOLAT[i-bk]);
         if(dsum>=limit)break;
        }

      int period=bk;
      for(int j=period-1;j<=InpMaxPeriod;j++)
        {
         if(S1[i-j]==EMPTY_VALUE)break;
         S1[i-j]=EMPTY_VALUE;
         R1[i-j]=EMPTY_VALUE;
         S2[i-j]=EMPTY_VALUE;
         R2[i-j]=EMPTY_VALUE;
        }
      //---
      double alpha=1,beta=0;
      calc_trend(alpha,beta,close,period,i);
      int from_x=i-(period-1);
      double y0 = beta-(period-1)*alpha;


      //---
      double upper[][2];
      double lower[][2];
      int up_n=0;
      int lo_n=0;
      copy_outside(upper,up_n,lower,lo_n,HI,LO,alpha,from_x,y0,i);

      //---
      double ra=0,rb=0;
      double r_ma=0;
      double r_max=0;
      calc_snr(ra,rb,r_ma,r_max,upper,up_n,i,1);
      //---
      double sa=0,sb=0;
      double s_ma=0;
      double s_max=0;
      calc_snr(sa,sb,s_ma,s_max,lower,lo_n,i,-1);

      double ry1,sy1,ry2,sy2;

      ry1=ra*from_x+rb+r_ma;
      sy1=sa*from_x+sb-s_ma;
      //---
      ry2=ra*from_x+rb+r_max;
      sy2=sa*from_x+sb-s_max;

      //---

      for(int j=from_x;j<=i;j++)
        {
         S1[j]=sy1;
         S2[j]=sy2;
         sy1+=sa;
         sy2+=sa;
         R1[j]=ry1;
         R2[j]=ry2;
         ry1+=ra;
         ry2+=ra;
        }


     }
//--- return value of prev_calculated for next call
   return(rates_total);
  }
//---
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
//|
//+------------------------------------------------------------------+
void calc_trend(double  &a,double  &b,const double &price[],const int cnt,const int to)
  {

   if(cnt==0)
     {
      a=EMPTY_VALUE;
      b=EMPTY_VALUE;
      return;
     }

   int from=to-(cnt-1);
   int x=0;
//--- 
   double sumy=0.0; double sumx=0.0;
   double sumxy=0.0; double sumx2=0.0;

   for(int k=from;k<=to;k++)
     {
      x++;
      sumx+=x;
      sumx2+=x*x;
      sumy+=price[k];
      sumxy+=x*price[k];
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
      b+=a*cnt;
     }
  }
//+------------------------------------------------------------------+
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
