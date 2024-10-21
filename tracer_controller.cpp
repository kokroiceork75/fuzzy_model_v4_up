// ~~~~~~~~~~~~~~~~~ 多目標改成單目標 cost function 合併 ~~~~~~~~~~~~~~~~~~~~
// modify this first
#include <assert.h>  
#include <math.h>
#include <stdio.h>  
#include <stdlib.h>   
#include <webots/lidar.h> 
#include <webots/motor.h>  
#include <webots/robot.h>   
#include <webots/supervisor.h>   
#include <iostream>
#include <cmath>
#include "time.h"   

#define _in_varl  2   
#define _out_varl 3
#define _conq_varl 9
#define ngen 600   /* 最大的世代數   */ 
#define popsize 40  /* 解的個數 */    ///202012
#define new_row 40  /* 新產生的解的個數  */	
#define _rule_number 10   /* 規則數 */ 
#define quarter_pop 12  /* 4分之1個解的數量，用在ACO解的更新  */	
#define  _max_rule  30    /* 開記憶體空間用的 */
#define nobj 1             /*  目標個數，目前為2表示是2個目標  */
#define   _pcross   0.9   /* 交配機率 */
#define   _pmutat   0.1  /* 突變機率 */
#define eta_c 10      /*  交配的變數  */
#define eta_m 10      /*  突變的變數  */
#define   garma02    0.6   //  調整可解釋性的合併情況，調整越大，合併條件越嚴格    
#define   garma03    0.75   //  調整可解釋性的costfunction的第一限制裡面的set之間的距離，調整越大，set會分的越散
#define   _conq_max  1.0       /*  後件部速度項最大值  */ // wei need change 範圍要乘以倍率 例如50
#define   _conq_min  -1.0      //  後件部速度項最小值
#define   _conq_const_max  5.0 // 後件部常數項最小值
#define   _conq_const_min  -5.0 // 後件部常數項最小值

#define   _range_max  0.2     /*  寬度最大值  */
#define   _range_min  0.05     /*  寬度最小值  */ 
#define   _center_max  1.0   /*  中心點最大值  */ // wei need change 倍率
#define   _center_min  -1.0    /*  中心點最小值   */ 
#define left_wheel_speed 50.0  /* 左輪輪速 */
#define right_wheel_speed  50.0   /* 右輪輪速  */
#define stop_max_dis 1.5    /* 停止的距離  */
#define stop_dis_narrow 0.5  /* ? */
#define dwall   0.535   /*  希望離牆的值  */ ///0.285(車體邊緣) + 0.35(離牆距離)
#define max_step 50      /* 最大步數  */  
#define max_tune_step 30   /*  重置點的步數  */
#define sensor_range 1 /* 感測的範圍:1m */
#define _input_scale (1/sensor_range)  /* 給感測器讀值正規化的參數  */
#define getArRadius 0   /*  車體半徑  */
#define max_step_to14 30
#define randomize() srand((unsigned) time(NULL)) /*  以時間為基底產生隨機值的初始化，不要動這行  */
#define random(x) rand() %x   /*  隨機產生亂數，不要動這行  */
#define TIME_STEP 32 /*  延遲時間 */
#define INF 1.0e14     /*  無限大  */
#define layer 8
#define angleeee 361

using namespace std ; 

double in[_in_varl+1] ; //感測器正規化後的輸入
double out_y[_out_varl+1] ; //馬達輸出 //0~5
double robot_x, robot_y, robot_fi;
double prev_robot_x, prev_robot_y, prev_robot_fi;
double delta_x, delta_y, delta_fi;
static int random_initialized = 0;
double sick_1[max_step+1],sick_2[max_step+1];
double delta_x_arr[max_step+1], delta_y_arr[max_step+1], delta_fi_arr[max_step+1];
double out_delta_x_arr[max_step+1], out_delta_y_arr[max_step+1], out_delta_fi_arr[max_step+1];

int step ;
int remained_step, current_stop_step ;
int c_gen;  //計數用的世代數,ngen是最大的世代數
int ber ; //代表第幾個解，此為一個替代的暫存變數
const int _mem_length = 2*_in_varl+ _conq_varl ; // 一個rule的長度 = 2*2 + 9 = 13
const int real_mem = _rule_number*(2*_in_varl+ _conq_varl) ;  // 整個解的個數 10*13=130
int _in_clu ; //rule_number的暫存變數	
double fuzzy_real[popsize+new_row+1][_mem_length*_max_rule+1] ; //  實際輸出規則的變數 0~80*0~300

//非主導排序相關
int *temp1=NULL ; //用在非主導排序的暫存變數
int front_size ; // 排除front1後，統計該回合產生多少個Front的轉存變數
int *temp_sort ; //用在非主導排序的暫存變數
int front_sort[2*popsize+1] ; //用在非主導排序的暫存變數
int sorted_crowding[2*popsize+1] ; //用在非主導排序的暫存變數 //0~80

//多樣性相關
double temp_crowd[2*popsize+1][nobj+2] ; //排列解的多樣性的暫存變數
int mixed_pop_sorted[2*popsize+1] ; //排列子代解的多樣性的暫存變數

double min_m[popsize+1][_in_varl+1],max_m[popsize+1][_in_varl+1];//最小中心點 & 最大中心點		
const double INITIAL[3] = { 0, 0.0891353, 0 };  //初始位置
float xx[max_step+1],yy[max_step+1]; //傳送座標資料時的暫存變數
int best_solution[nobj+1]={0};
double fuzzy_real_compare[_mem_length*_max_rule+1] ; //儲存"已排除重複規則"的暫存空間
int load_gen;  	
long double parent_rule_load[popsize+1][real_mem+1]={0.0};
long double parent_obj_load[popsize+1][nobj+1]={0.0};
int parent_rank_load[popsize+1]={0};

void robot_position(int);  //抓取機器人位置
void robot_position_reset(int);  //重置機器人位置
void evaluate_pop (int,int,int,int,int) ; // 計算pop:用來評分、計算cost value與限制population
void initialize_pop (); // 將解初始化
void update_newparent() ;    //更新親代
int tour_selection() ; //競爭法
void merge(); //混合
void update_pop(double , double , int) ;   //解的更新
void crossover(int , int , int) ; //交配
void mutation() ;    // 突變  
void check_and_combime_child_set(int);  //fuzzy set的合併
void child_pop_cluster(int,int,int); //子代群
void overlap_degree(int,int);  //重疊度的計算
void fuzzy_set_constraint(int,int);  //fuzzy set之間的距離限制
void non_domination(int,int) ;  //非主導排序:排列解的好壞
void sort_on_obj(int,int) ;  //按照 cost value 將解排序
void sort_on_crowd(int , int) ; //排列解的多樣性
void sort_on_rank_crowd(int) ; //做排列子代解的多樣性
void In_mem_fir(double *,int,int) ; //計算激發量
void Fir_str(double *, int,int) ; //Fuzzy系統:讀取後件部與計算激發量
void fuzzy(int,int,int,double,double);    // Fuzzy系統:解模糊
void fuzzy_in(int) ;  // Fuzzy系統:把感測器的讀值做正規化
void report_pop (FILE *fpt);   //回傳解
void report_best(int,int) ; //回傳最佳解
void report_pop_for_reading ( FILE *fpt) ;
void save_final_obj(int) ; //存最終fuzzy set解
void save_time(float total_time) ;  //儲存時間      
void save_final_pareto_front() ;  //存最後一代的 front 
void save_temp_obj(int) ; //存暫時fuzzy set解
void gauss_data_limit(int) ;  // 限制高斯函數資料範圍         
void get_sensor(int,float*) ; //抓感測器讀值的            
void clear_matrix(void) ; //清除矩陣資料
// void speed_limit(void) ;   // 速度限制器的原型宣告
void data_printf(int,int,int,int) ;  //印資料
double randn(double, double); //產生隨機亂數
void print_gen_obj(int);  //印出親代的cost
void execute_process(int);  //主執行序
void load_gen_data();
int end_iteration_and_save_gen(int);
void load_save_parent_rule_and_relative_data();

float get_random_velocity();
void calculate_robot_delta(int);
void cost_function(int);

class Population
{
friend class C_Rule ;
public :
 
double xreal[_max_rule*_mem_length];
double crowd_dist;
double objv[nobj+1] ;   //此為cost function的數量
double smoothly[2];
int sp[2*popsize+1] ;
int nsp ;
int sorted_front ;
double dist[nobj+2] ;
int sorted_crowd ;
int nnp ;
int np[2*popsize+1] ;
int front_rank ;
Population() ;
} ;

Population :: Population()
{
}

Population parent_pop[popsize+new_row+1] ; //親代空間 0~80
Population child_pop[popsize+new_row+1] ; // 子代空間 0~80
Population mixed_pop[2*popsize+1] ; // 混合與後續功能用來轉存的暫存變數 0~80
Population temp_pop[2*popsize+1] ; //非主導排序與後續功能用來轉存的暫存變數 0~80

class Front
{
public:
  int member[2*popsize+1] ; //0~80
  int number ;
} ;

Front _front[2*popsize+1] ;

inline double phi(double x,double m,double v)
{ return(  exp( - (x-m)*(x-m)/(v*v) )); }

class C_Rule
{
friend class Population ;
public :

double  in_mu ;
double  con[_conq_varl+1] ;
void    In_mem_fir(double *,int,int) ;
friend  void Fir_str(double *, int,int) ;
// friend  void fuzzy(double *,int,int); // wei need change,seems not use? change 註釋
} ;

C_Rule  _Rule[_max_rule+1] ;

float min_variable(int i) // 設置前件部、後件部參數的最小值
{
  float ge;
  int j;
  int n=2*_in_varl + _conq_varl;   // n = 13
  j = (i - 1) % n + 1; // 將 i 映射到 [1, 13] 範圍
  
  if (j == 1 || j == 3) {
      ge = _center_min;
  } else if (j == 2 || j == 4) {
      ge = _range_min;
  } else if (j == 5 || j == 6 || j == 8 || j == 9 || j == 11 || j == 12) {
      ge = _conq_min;
  } else { // 7,10,13 
      ge = _conq_const_min;
  }
  return ge;
}
float max_variable(int i) // 設置前件部、後件部參數的最小值
{
  float ge;
  int j;
  int n=2*_in_varl + _conq_varl;   // n = 13
  j = (i - 1) % n + 1; // 將 i 映射到 [1, 13] 範圍
  
  if (j == 1 || j == 3) {
      ge = _center_max;
  } else if (j == 2 || j == 4) {
      ge = _range_max;
  } else if (j == 5 || j == 6 || j == 8 || j == 9 || j == 11 || j == 12) {
      ge = _conq_max;
  } else { // 7,10,13 
      ge = _conq_const_max;
  }
  return ge;
}


int main(int argc, char **argv) 
{
  load_gen_data();
  if(load_gen <= ngen)
  {
    if(load_gen==1)
    {
      execute_process(load_gen);  //主執行序
      load_gen = end_iteration_and_save_gen(load_gen);
      printf("final=%d\n",load_gen);
      wb_supervisor_world_reload();
      wb_robot_step(TIME_STEP);    
    }
    else
    {
      load_save_parent_rule_and_relative_data();    //親代參數
      execute_process(load_gen);   //主執行序
    
      load_gen = end_iteration_and_save_gen(load_gen);

      printf("final=%d\n",load_gen);

      wb_supervisor_world_reload();
      wb_robot_step(TIME_STEP); 

    }

  }
  wb_robot_cleanup();
  return 0;
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void execute_process(int c_gen)
{
  int i = 0, j, kk, sw, jj, jjj;
  int iterate_counter;
  int front_one_counter=0; //計算 front one 的個數的計數器
  FILE *fpt4,*fpt6;  //提供存檔的指標變數
  randomize() ;  // 以時間為基底產生隨機值的初始化 // wei need change
  wb_robot_init();
  WbDeviceTag left_wheel = wb_robot_get_device("left_motor");
  WbDeviceTag right_wheel = wb_robot_get_device("right_motor");

  wb_motor_set_position(left_wheel, INFINITY);
  wb_motor_set_position(right_wheel, INFINITY);
  
  double left_speed = 0.0, right_speed = 0.0;
  wb_motor_set_velocity(left_wheel, left_speed);
  wb_motor_set_velocity(right_wheel, right_speed);

  _in_clu = 0 ; //rule_number的暫存變數
  if(c_gen == 1)
  {
    initialize_pop(); //將解初始化
    for(int ss=1; ss<=popsize; ss++ )
    {
      for(int ssss=1; ssss<=_in_varl; ssss++)
      {
        min_m[ss][ssss]=1;  
        max_m[ss][ssss]=0;  
      }
    }

    for(i=1; i<=popsize; i++ )
    {
      for(jj=1; jj<=_rule_number; jj++)
      {
          for(jjj=1; jjj<=_in_varl; jjj++)
          {       
            // 前件部轉存 
            fuzzy_real[i][(jj-1)*_mem_length + jjj*2 - 1] = parent_pop[i].xreal[(jj-1)*_mem_length+jjj*2 - 1] ; // 中心點
            fuzzy_real[i][(jj-1)*_mem_length + jjj*2] = parent_pop[i].xreal[(jj-1)*_mem_length + jjj*2] ; // 寬度
            if(fuzzy_real[i][(jj-1)*_mem_length+jjj*2-1] < min_m[i][jjj])
            {
              min_m[i][jjj]=fuzzy_real[i][(jj-1)*_mem_length+jjj*2-1];  //如果比最大的最小值(min_m = 1)還小，就儲存最小中心點
            }
            
            if(fuzzy_real[i][(jj-1)*_mem_length+jjj*2-1] > max_m[i][jjj])
            {
              max_m[i][jjj]=fuzzy_real[i][(jj-1)*_mem_length+jjj*2-1];  //如果比最小的最大值(max_m = 1)還大，就儲存最大中心點，
            }
              
          }
          // fuzzy_real[i][jj*_mem_length-_out_varl+1] = parent_pop[i].xreal[jj*_mem_length-_out_varl+1] * left_wheel_speed ; // 左輪
          // fuzzy_real[i][jj*_mem_length-_out_varl+2] =  parent_pop[i].xreal[jj*_mem_length-_out_varl+2] * right_wheel_speed ;  // 右輪
          // 後件部轉存
          for(int k=1; k<=9; k++)
          {
            fuzzy_real[i][jj*_mem_length-_conq_varl + k] =  parent_pop[i].xreal[jj*_mem_length - _conq_varl + k]; 
          }
            
          
      }   
    }

    iterate_counter = 0;
    ber = 0;
    while (ber<popsize) 
    {
      // int lser=0; ///202011 為了防止讀到前一個控制器的雷射值
      ber=ber+1;
      parent_pop[ber].objv[1]=0; // 將親代的第一個 cost value歸零。如果沒有5位置訓練要先歸零(有的話則不用)
      parent_pop[ber].objv[2]=0; // 將親代的第二個 cost value歸零。如果沒有5位置訓練要先歸零(有的話則不用)
      // fail_step_record[ber]=0;  //如果沒有5位置訓練要先歸零(有的話則不用)
      // total_sucess_step[ber]=0; //如果沒有5位置訓練要先歸零(有的話則不用)
      // total_fail_step=0;  //失敗的步數歸零  
      printf("ber = %d\n",ber);
      // wei need change 
      // int step=0; // 原本是沒註解，為了fuzzy()改全局
      step = 0; // wei need change first, not useful
      robot_position(step);
      remained_step = 0 ; //殘留步數
      current_stop_step = max_step ;   //max_tune_step:重置點的步數 
      //停止的步數=重置點的步數
      _in_clu = _rule_number;  //rule_number的暫存變數
      int far_reset_counter = 0; //距離過遠而重置的計數器
      while (step < max_step) 
      { 
        step=step+1;  
        printf("親代執行區!!!!!!!\n");
        wb_robot_step(100);
        robot_position(step);
        // wei need change 這裡給予隨機速度
        sick_1[step] = get_random_velocity();
        sick_2[step] = get_random_velocity();
        printf("速度 sick_1: %lf\n",sick_1[step]); // wei debug
        printf("速度 sick_2: %lf\n",sick_2[step]);
        double vl = sick_1[step];
        double vr = sick_2[step];
        wb_motor_set_velocity(left_wheel, vl);
        wb_motor_set_velocity(right_wheel, vr);
        fuzzy_in(step) ;  //輸入是雷射數值，只是調整fuzzy的input比例(做正規化)			
        Fir_str(in , _rule_number , ber) ; // 讀取後件部和計算每條rule激發量
        fuzzy(_rule_number, ber, step, vl, vr); // 解模糊產生out_y
        // speed_limit() ;
        data_printf(step, c_gen, ber, _in_clu) ; //印資料	
        // wb_motor_set_velocity(left_wheel , out_y[1]);
        // wb_motor_set_velocity(right_wheel, out_y[2]);
        // 再次記錄位置 robot_position
        wb_robot_step(100);
        robot_position(step);
        // 計算這50毫秒的姿態變化量
        calculate_robot_delta(step);

        // wei need change 需要計算sigma cost value 
  
        // if ( far_reset_counter>20 || margin_reset==1 || no_move_reset_counter == 200 || (out_y[1] <= 0.0 && out_y[2] <= 0.0))
        // { ... }			

        if(step==max_step)
        {
          wb_motor_set_velocity(left_wheel , 0); 
          wb_motor_set_velocity(right_wheel, 0);
          printf(" 停止，收集到一定數據 \n") ;        
          break;
        }
        fflush(stdout);
      }
      cost_function(c_gen);
      // printf("~~~~~~~~~~~~~~~~~~cost function 1: %lf, cost function 2: %lf~~~~~~~~~~~~~~~~~~\n", f1, f2);
      // evaluate_pop (c_gen, ber, _in_clu, remained_step, current_stop_step) ;  //計算costvalue的前置作業 	
      clear_matrix() ; // 清除陣列資料

      // parent_pop[ber].objv[1] = (parent_pop[ber].objv[1] / total_sucess_step[ber]) + fail_step_record[ber];   
      // parent_pop[ber].objv[2] = ( 1/(parent_pop[ber].objv[2] / total_sucess_step[ber]) )+ fail_step_record[ber];  
      print_gen_obj(c_gen);  //印出親代的cost
      robot_position_reset(0);  

      //overlap_degree(ber, _in_clu);  //重疊度的計算
      //fuzzy_set_constraint(ber, _in_clu);  //fuzzy set之間的距離限制 
      if(ber==popsize)
      {
        break;
      }
     
    }

    sw = 1 ; // 使親代進入非主導排序的條件 
    non_domination(sw , c_gen) ; //非主導排序
    sort_on_crowd(sw , c_gen) ; //排列解的多樣性
    printf("\n gen = %d\n", c_gen); //印出目前是第幾代
    report_best(c_gen,sw) ; //回傳最佳解

    FILE *fnoise1;  
    fnoise1 = fopen("parent_rank.txt","w"); 	
    for (i=1; i<=popsize; i++)
    {
      fprintf(fnoise1,"%d\t", parent_pop[i].front_rank); //排名
      fprintf(fnoise1, "\n") ;
    } 
    fclose(fnoise1);

    FILE *fnoise2;
    fnoise2 = fopen("parent_obj.txt","w"); 	
    for (i=1; i<=popsize; i++)
    {
      for(j=1; j<=nobj;j++)
      {
        fprintf(fnoise2,"%f\t", parent_pop[i].objv[j]) ;  //cost value
      }
      fprintf(fnoise2, "\n") ;
    }
    fclose(fnoise2);


    FILE *fnoise3; 
    fnoise3 = fopen("parent_rule.txt","w"); 	
    for (i=1; i<=popsize; i++)
    {   
      for(jj=1; jj<=_rule_number; jj++)
      {
        for(jjj=1; jjj<=_in_varl; jjj++)
        {
          fprintf(fnoise3,"%.10f\t", parent_pop[i].xreal[(jj-1)*_mem_length+jjj*2-1] ); // 中心點
          fprintf(fnoise3,"%.10f\t", parent_pop[i].xreal[(jj-1)*_mem_length+jjj*2] ) ; // 寬度
        }
        // fprintf(fnoise3,"%.10f   \t", parent_pop[i].xreal[jj*_mem_length-_out_varl+1] ) ; // 左輪
        // fprintf(fnoise3,"%.10f   \t", parent_pop[i].xreal[jj*_mem_length-_out_varl+2] ) ;  // 右輪  	
        for(int k =1; k <= _conq_varl; k++)		
        {
          fprintf(fnoise3,"%.10f\t", parent_pop[i].xreal[jj*_mem_length - _conq_varl + k] ) ;  // 後件部 wei debug	
        }
      }    
      fprintf(fnoise3, "\n") ;
    }
    fclose(fnoise3);

  // c_gen = 1 結束
  }
  else if(c_gen>1)
  {
    printf("I'm the %d gen!!!!!!!\n",c_gen);
    for (i=1; i<=popsize; i++)
    {
      for (j=1; j<=_max_rule*_mem_length; j++)
      {
         fuzzy_real[i][j] = 0. ;
      }
    }
    double _cg = c_gen ;
    double _mg = ngen ;
    int t_loop = new_row; //做為轉存變數 t_loop = 40
    kk=1;
    update_pop(_cg , _mg , kk) ; // 解的更新
    
    for (kk=1; kk<=new_row; kk++) 
    {
      gauss_data_limit(kk); // 限制高斯函數資料範圍 
    }

	/////////////////  最小中心點 & 最大中心點 的初始化  ///////////////{

	for(int ss=1; ss<=new_row; ss++ )
	{
		for(int ssss=1; ssss<=_in_varl; ssss++)
		{
			min_m[ss][ssss]=1;  //為了找最小值，所以初始的最小值令為1，即為"最大"的最小值。
			max_m[ss][ssss]=0;  //為了找最大值，所以初始的最大值令為0，即為"最小"的最大值。

			//驗證用
			//printf("min_m[%d][%d]:%f\n",ss,ssss,min_m[ss][ssss]); 
		}
	}
  iterate_counter=0;
  ber=0;
  while (ber<new_row) 
  {
    // int lser=0; 
    // while (lser<2)
    // {
    //   lser=lser+1;
    //   printf("雷射讀值!!!!!!!\n");
    //   wb_robot_step(10);
    //   const float *lidar_temp_value = wb_lidar_get_layer_range_image(lidar, layer);
    //   float *lidar_value = new float[361];
    //   for (int i = 0; i <= horizontalResolution+1; ++i) 
    //   {
    //     int angle = i * 360.0 / (horizontalResolution+1);
    //     lidar_value[angle] = lidar_temp_value[i];
    //   }
    //   get_sensor(step,lidar_value);
    //   wb_motor_set_velocity(left_wheel , 0);
    //   wb_motor_set_velocity(right_wheel, 0);
    //   wb_robot_step(100); 
    //   if(lser==2) { break; }
    // }
    ber=ber+1;
    child_pop[ber].objv[1]=0; 
    child_pop[ber].objv[2]=0; 
    // fail_step_record[ber]=0;  
    // total_sucess_step[ber]=0; 
    // total_fail_step=0;  
    printf("ber = %d\n",ber);
    int step=0;
  	robot_position(step);
    remained_step = 0 ; //殘留步數
    current_stop_step = max_step ;   //max_tune_step:重置點的步數 
    _in_clu = _rule_number;  //rule_number的暫存變數
    int far_reset_counter = 0; //距離過遠而重置的計數器
    double narrow_down = 0.; //?
    while (step<max_step) 
    { 
      step=step+1;  	
      printf("子代執行區!!!!!!!\n");
      wb_robot_step(100);
      robot_position(step);
      // wei need change 這裡給予隨機速度
      sick_1[step] = get_random_velocity();
      sick_2[step] = get_random_velocity();
      double vl = sick_1[step];
      double vr = sick_2[step];
      wb_motor_set_velocity(left_wheel, vl);
      wb_motor_set_velocity(right_wheel, vr);
      fuzzy_in(step) ;  
		  child_pop_cluster(ber,step,c_gen); // 跟親代不一樣
      Fir_str(in , _rule_number , ber) ; //讀取後件部計算每條rule機發量
      fuzzy(_rule_number , ber, step, vl, vr); 
      // speed_limit() ; // 速度限制器
      data_printf(step, c_gen, ber, _in_clu) ; //印資料	
      // wb_motor_set_velocity(left_wheel , out_y[1]);
      // wb_motor_set_velocity(right_wheel, out_y[2]);
      wb_robot_step(100);
      robot_position(step);
      calculate_robot_delta(step);
      // narrow_down 似乎完全沒用到
      // if(step >= 0.25*max_step && step <= 0.75*max_step)  
      // {
      //    narrow_down = narrow_down + (stop_max_dis-stop_dis_narrow)/(0.5*max_step);
      //    //printf("stop_max - narrow_down = %.8f\n",stop_max_dis-narrow_down);
      // }    
      // if(sick_right[step] > (stop_max_dis - narrow_down))
      // {
      //     far_reset_counter++;
      // }
      // else
      // {
      //     far_reset_counter = 0; //離牆太遠的計數器歸零
      // }
      
      // if ( far_reset_counter>20 || margin_reset==1 || no_move_reset_counter == 200 ||(out_y[1] <= 0.0 && out_y[2] <= 0.0))
      // {......     }			
		  if(step==max_step)
		  {
			  wb_motor_set_velocity(left_wheel , 0); 
        wb_motor_set_velocity(right_wheel, 0);
        printf(" 停止，收集到一定數據 \n") ;        
        break;
		  }
      fflush(stdout);  
    }
    printf("**************************%d  %d\n", c_gen, _cg);
    cost_function(c_gen);
    // evaluate_pop (c_gen, ber, _in_clu, remained_step, current_stop_step) ;  //計算costvalue的前置作業 	
    clear_matrix() ; // 清除陣列資料
    // child_pop[ber].objv[1] = (child_pop[ber].objv[1] / total_sucess_step[ber]) + fail_step_record[ber];
    // child_pop[ber].objv[2] = (1/(child_pop[ber].objv[2] / total_sucess_step[ber])) + fail_step_record[ber];
    print_gen_obj(c_gen);  //印出親代的cost
    robot_position_reset(0);
	  // printf("total_sucess_step[%d]= %f\n",ber,total_sucess_step[ber]);
	  // printf("fail_step_record[%d]= %f\n",ber,fail_step_record[ber]);
	  // child_pop[ber].objv[1] = child_pop[ber].objv[1] ;// + 0.5*(child_pop[ber].smoothly[1]);
	  // child_pop[ber].objv[2] = child_pop[ber].objv[2] ;// + 0.25*(child_pop[ber].smoothly[1]);
	  // overlap_degree(ber, _in_clu);  //重疊度的計算
	  // fuzzy_set_constraint(ber, _in_clu);  //fuzzy set之間的距離限制
    if(ber==new_row)
    {
      break;
    }   
  }    
  merge();  // 親袋子代一起放入mix_pop
  sw = 3 ;// 使子代進入非主導排序的條件  
  non_domination(sw , c_gen) ; //非主導排序
  sort_on_crowd(sw , c_gen) ;  //排列解的多樣性
  update_newparent() ; //使用參數歸零
  report_best(c_gen,sw) ; //回傳最佳解
 
  FILE *fnoise1;  
  fnoise1 = fopen("parent_rank.txt","w"); 	
  for (i=1; i<=popsize; i++)
  {
    fprintf(fnoise1,"%d\t", parent_pop[i].front_rank); //排名
    fprintf(fnoise1, "\n") ;
  }
  fclose(fnoise1);

FILE *fnoise2;
  fnoise2 = fopen("parent_obj.txt","w"); 	
  for (i=1; i<=popsize; i++)
  {
    for(j=1; j<=nobj;j++)
    {
      fprintf(fnoise2,"%f\t", parent_pop[i].objv[j]) ;  //cost value
    }
     fprintf(fnoise2, "\n") ;
  }
  fclose(fnoise2);


FILE *fnoise3; 
   fnoise3 = fopen("parent_rule.txt","w"); 	
  for (i=1; i<=popsize; i++)
  {   
    for(jj=1; jj<=_rule_number; jj++)
    {
      for(jjj=1; jjj<=_in_varl; jjj++)
      {
        fprintf(fnoise3,"%.10f\t", parent_pop[i].xreal[(jj-1)*_mem_length+jjj*2-1] ); // 中心點
        fprintf(fnoise3,"%.10f\t", parent_pop[i].xreal[(jj-1)*_mem_length+jjj*2] ) ; // 寬度
      }
    //  fprintf(fnoise3,"%.10f   \t",  parent_pop[i].xreal[jj*_mem_length-_out_varl+1] ) ; // 左輪
    //  fprintf(fnoise3,"%.10f   \t",  parent_pop[i].xreal[jj*_mem_length-_out_varl+2]  ) ;  // 右輪  
    for(int k =1; k <= _conq_varl; k++)		
    {
      fprintf(fnoise3,"%.10f\t", parent_pop[i].xreal[jj*_mem_length - _conq_varl + k] ) ;  // 後件部  wei debug
    }			
    }    
    fprintf(fnoise3, "\n") ;
  }
  fclose(fnoise3);
  
  		for(int j=1; j<=nobj;j++)
		{
			save_temp_obj(j) ; //存最終fuzzy set解
		}	
   // Save all pop
  fpt4 = fopen("..\\..\\generated_file\\all_pop.out","a");
  fprintf(fpt4,"# gen = %d\n ", c_gen);  
  report_pop(fpt4);  // all_pop.out
  fflush(stdout);  // fflush() 就是在「強迫」他將 buffer 中的東西輸出到螢幕.
  fflush(fpt4);
  fclose(fpt4);
  fflush(stdout);

  // c_gen > 1結束
 } 
  if(c_gen==ngen)  
  {  
    printf("\n Generations finished, now reporting solutions\n");
    save_final_pareto_front() ; //存最後一代的 front 
    for(int j=1; j<=nobj;j++)
    {
      save_final_obj(j) ; //存最終fuzzy set解
    }	
  }
  // execute_process() 結束
}

void initialize_pop() // wei初始化所有規則參數
{
  int i, j, jj, jjj;
  for (i=1; i<=popsize; i++)
  {
    for (j=1; j<=_max_rule*_mem_length; j++)
    {
      fuzzy_real[i][j] = 0. ;
    }
  }

  for (i=1; i<=popsize; i++ )
  {
    for (jj=1; jj<=_rule_number; jj++)
    {			
      for(jjj=1; jjj<=_in_varl; jjj++)
      {
        parent_pop[i].xreal[(jj-1)*_mem_length + jjj*2 - 1] = 1.0*(double)rand()/RAND_MAX  ;  // 中心點
        parent_pop[i].xreal[(jj-1)*_mem_length + jjj*2] =  _range_min + (fabs(_range_max-_range_min))*(double)rand()/RAND_MAX ; // 寬度[_range_min, _range_max]    
      }		
      // 後件部隨機初始化		
      parent_pop[i].xreal[jj*_mem_length - _conq_varl + 1] = 50. * ( 2.0*(double)rand()/RAND_MAX -1.0 ) ; // wei need change 
      parent_pop[i].xreal[jj*_mem_length - _conq_varl + 2] = 50. * ( 2.0*(double)rand()/RAND_MAX -1.0 ) ;  // [-50~50]隨機浮點數
      parent_pop[i].xreal[jj*_mem_length - _conq_varl + 4] = 50. * ( 2.0*(double)rand()/RAND_MAX -1.0 ) ; // conq v
      parent_pop[i].xreal[jj*_mem_length - _conq_varl + 5] = 50. * ( 2.0*(double)rand()/RAND_MAX -1.0 ) ; 
      parent_pop[i].xreal[jj*_mem_length - _conq_varl + 7] = 50. * ( 2.0*(double)rand()/RAND_MAX -1.0 ) ; 
      parent_pop[i].xreal[jj*_mem_length - _conq_varl + 8] = 50. * ( 2.0*(double)rand()/RAND_MAX -1.0 ) ; 
      parent_pop[i].xreal[jj*_mem_length - _conq_varl + 3] = 5. * ( 2.0*(double)rand()/RAND_MAX -1.0 ) ; // wei need change
      parent_pop[i].xreal[jj*_mem_length - _conq_varl + 6] = 5. * ( 2.0*(double)rand()/RAND_MAX -1.0 ) ; // [-5~5]隨機浮點數
      parent_pop[i].xreal[jj*_mem_length - _conq_varl + 9] = 5. * ( 2.0*(double)rand()/RAND_MAX -1.0 ) ; // conq const
    }
  }

}

void robot_position(int k) // 抓真實機器人位置和角度
{
  prev_robot_x = robot_x;
  prev_robot_y = robot_y;
  prev_robot_fi = robot_fi;

  WbNodeRef robot_node = wb_supervisor_node_get_from_def("Tracer");
  WbFieldRef trans_field = wb_supervisor_node_get_field(robot_node, "translation");
  const double *trans = wb_supervisor_field_get_sf_vec3f(trans_field);
  WbFieldRef rot_field = wb_supervisor_node_get_field(robot_node, "rotation");
  const double *rot = wb_supervisor_field_get_sf_rotation(rot_field);
  double yaw;
  // yaw 會存基於世界座標的車頭弧度 wei need change
  if (rot[1] > 0.99) {
      // 當旋轉軸幾乎與 y 軸平行時
      yaw = 2 * atan2(rot[0], rot[3]);
  } else if (rot[1] < -0.99) {
      // 當旋轉軸幾乎與 -y 軸平行時
      yaw = -2 * atan2(rot[0], rot[3]);
  } else {
      // 一般情況
      yaw = atan2(2 * (rot[3] * rot[1] - rot[0] * rot[2]), 
                  1 - 2 * (rot[1] * rot[1] + rot[2] * rot[2]));
  }

  // 調整弧度範圍到 [0, 2π)
  if (yaw < 0) {
      yaw += 2 * M_PI;
  }

  robot_x = trans[0];
  robot_y = trans[2];
  robot_fi = yaw;
 fflush(stdout);
}

void robot_position_reset(int state)
{
  WbNodeRef robot_node = wb_supervisor_node_get_from_def("Tracer");
  WbFieldRef trans_field = wb_supervisor_node_get_field(robot_node, "translation");
  WbFieldRef  rot_field  =  wb_supervisor_node_get_field(robot_node,"rotation");
  printf(" 重置! \n") ;
  if(state==0) // 可以設置不同重置點
   {  
     const double INITIAL[3] = { 0, 0.0891353, 0 };
     const  double  INITIAL_ROT[4]  =  { 0,  1, 0, 1.5708 };
     wb_supervisor_field_set_sf_vec3f(trans_field, INITIAL);   
     wb_supervisor_field_set_sf_rotation(rot_field,  INITIAL_ROT);
     wb_supervisor_simulation_reset_physics();
   }
   fflush(stdout);
}


float get_random_velocity()
{
    if (!random_initialized) {
        srand(time(NULL));
        random_initialized = 1;
    }
    
    // 生成 0 到 1 之間的隨機浮點數
    float random_float = (float)rand() / RAND_MAX;
    
    // 將範圍調整到 -50 到 50
    return random_float * 100.0f - 50.0f;
}

void fuzzy_in(int jj) // 正規化輸入
{
  in[1] = sick_1[jj] / 50.0; // wei need change, normalization
  in[2] = sick_2[jj] / 50.0;
}

void Fir_str(double *in, int _rule, int ber)
{
  int j,k ;

  for(k=1; k<=_rule; k++)  //_rule值=_rule_number=10
  {
    for (j=1; j<=_conq_varl; j++)
    {
      _Rule[k].con[j] = fuzzy_real[ber][k*_mem_length - _conq_varl + j] ; // 將後件部值做"轉存"
    } 	             
  _Rule[k].In_mem_fir(in,k,ber) ; // 計算第k條規則的激發量  
  }
}

void C_Rule::In_mem_fir(double *in, int _rule, int ber)
{
  // _rule代表第k條規則
  int i ;
  in_mu =  1. ;  // "1." 代表1的float值
  for(i=1;i<= _in_varl;i++)
  {
    if(in[i] < min_m[ber][i] || in[i] > max_m[ber][i])  
    //如果輸入的感測器正規化小於"最大的中心點" 或是 輸入的感測器正規化大於"最小的中心點"
    {
      in_mu = in_mu * 1;    //保持原樣，所以乘以1
    }
    else
    { 
      in_mu = in_mu * phi(in[i] , fuzzy_real[ber][(_rule-1)*_mem_length + i*2 - 1] , fuzzy_real[ber][(_rule-1)*_mem_length + i*2] ) ;
    }
  }
}


void fuzzy(int _rule, int ber, int step, double vl, double vr) 
{
    int i, j;
    double den[_out_varl+1],   //分母項 0~4
           num[_out_varl+1];   //分子項 0~4

    for (j = 1; j <= _out_varl; j++)
    {
        den[j] = 0.0;  //初始化歸零
        num[j] = 0.0;  //初始化歸零

        for (i = 1; i <= _rule; i++)
        {
            // 後件部使用矩陣運算，並從 _Rule[i].con 中提取對應的參數
            // [delta_x, delta_y, delta_fi] = [con[1] con[2] con[3]; con[4] con[5] con[6]; con[7] con[8] con[9]] * [vl vr 1]
            
            // 計算 delta_x
            num[1] += _Rule[i].in_mu * (_Rule[i].con[1] * vl + _Rule[i].con[2] * vr + _Rule[i].con[3]); 
            // 計算 delta_y
            num[2] += _Rule[i].in_mu * (_Rule[i].con[4] * vl + _Rule[i].con[5] * vr + _Rule[i].con[6]); 
            // 計算 delta_fi
            num[3] += _Rule[i].in_mu * (_Rule[i].con[7] * vl + _Rule[i].con[8] * vr + _Rule[i].con[9]); 

            // 計算分母部分，即各規則的激發量總和
            den[j] += _Rule[i].in_mu;
        }

        if (fabs(den[j]) < 1E-8)
            out_y[j] = 0;     //如果den很小時，直接設為0
        else
            out_y[j] = num[j] / den[j];  //加權平均法
        
        out_delta_x_arr[step] = out_y[1];
        out_delta_y_arr[step] = out_y[2];
        out_delta_fi_arr[step] = out_y[3];
        
    }
}

void data_printf(int step, int c_gen, int ber, int _in_clu)
{
	printf(" ~~~~~~~~~~~~~~~~~~~~~~~~~~ \n" ) ;
	printf("iteration == %d\n", c_gen) ;
	printf("ber == %d\n", ber) ;
	printf("time step == %d\n", step) ;
	printf("Sick_1 == %f\n", sick_1[step] ) ;
	printf("Sick_2 == %f\n", sick_2[step] ) ;
	printf("變化量x == %f\n", out_y[1]) ;
	printf("變化量y == %f\n", out_y[2]) ; 
	printf("變化量fi == %f\n", out_y[3]) ; 
	// printf("parent_real == %f\n", pa) ; 

	fflush(stdout);	
}

void calculate_robot_delta(int step) {
    delta_x = prev_robot_x - robot_x;
    delta_y = prev_robot_y - robot_y;
    delta_fi = robot_fi - prev_robot_fi;  
    delta_fi = fmod(delta_fi + M_PI, 2 * M_PI) - M_PI; // 調整 delta_fi 使其落在 -π 到 π 的範圍內
    delta_x_arr[step] = delta_x;
    delta_y_arr[step] = delta_y;
    delta_fi_arr[step] = delta_fi;
    

}
void evaluate_pop (int c_gen, int ber, int _rule, int remained_step , int current_stop_step)
{
  // rule = 10

}

void cost_function(int c_gen)
{
  // wei need change, c_gen and ber are global variables
  double f1 = 0.0, f2 = 0.0;
  double error_x = 0.0, error_y = 0.0, error_fi = 0.0;
  for(int i = 1; i <= max_step; i++)
  {
    error_x += pow(delta_x_arr[i] - out_delta_x_arr[i], 2);
    error_y += pow(delta_y_arr[i] - out_delta_y_arr[i], 2);
    error_fi += pow(delta_fi_arr[i] - out_delta_fi_arr[i], 2);

  }
  f1 = error_x + error_y + error_fi;
  // f2 = error_fi;

  if(c_gen > 1)
  {
    child_pop[ber].objv[1] = f1 + child_pop[ber].objv[1];
    child_pop[ber].objv[2] = f2 + child_pop[ber].objv[2];
    printf("OOOOOOOOOOOOOOOOOOOOOOO\n");
    printf("1111111After update: child_pop[ber].objv[1] = %f, child_pop[ber].objv[2] = %f\n", child_pop[ber].objv[1], child_pop[ber].objv[2]);
  }
  else // c_gen = 1
  {
    parent_pop[ber].objv[1] = f1 + parent_pop[ber].objv[1];
    parent_pop[ber].objv[2] = f2 + parent_pop[ber].objv[2];
    printf("OOOOOOOOOOOOOOOOOOOOOOO\n");
    printf("OOOOOOOOOOOOOOOOOOOOOOO\n");
    printf("OOOOOOOOOOOOOOOOOOOOOOO\n");
  }
  printf("Generation: %d, ber: %d, f1: %f, f2: %f\n", c_gen, ber, f1, f2); // wei debug
  printf("Updating objv for gen: %d, ber: %d, f1: %f, f2: %f\n", c_gen, ber, f1, f2);
  printf("Updated child_pop objv[1]: %f, objv[2]: %f\n", child_pop[ber].objv[1], child_pop[ber].objv[2]);

}

void clear_matrix() // 清除一些全局變量
{
  for(int i = 1; i <= max_step; i++)
  {
    sick_1[i] = 0.0;
    sick_2[i] = 0.0;
    delta_x_arr[i] = 0.0;
    delta_y_arr[i] = 0.0;
    delta_fi_arr[i] = 0.0;
    out_delta_x_arr[i] = 0.0;
    out_delta_y_arr[i] = 0.0;
    out_delta_fi_arr[i] = 0.0;

  }
}
void print_gen_obj(int gen)
{
   if(gen==1)  
   {
     printf("cost value 1 = %f\n",parent_pop[ber].objv[1]);    
     printf("cost value 2 = %f\n",parent_pop[ber].objv[2]); 
     printf("ccccccccccccccccccccccccccc_gen：　%d\n \n \n \n", c_gen);
   }    
   else if(gen>1) 
   {
      printf("cost value 1 = %f\n",child_pop[ber].objv[1]);
      printf("cost value 2 = %f\n",child_pop[ber].objv[2]);  
     printf("ccccccccccccccccccccccc_gen:　%d\n\n\n\n", c_gen);
   }
 fflush(stdout);
}

// wei need change the following all codes
void non_domination(int sw , int gg)
{
  // gg = c_gen
  int i,j,k,w,m ;
  double _fmin=0. , _fmax=0., next_obj=0. , pre_obj=0. ;
  double _dist=0.;
  int less, equal, more ;
  int front_rank ; 
  int pop_size=0 ;
  if (sw == 1)
  {
    pop_size = popsize ; //轉存解的數量
    
    for(i=1; i<=pop_size; i++)
    {
      for (j=1; j<=nobj; j++)
      {
        temp_pop[i].objv[j] = parent_pop[i].objv[j] ; //cost value的轉存
      }
    }
  }

  else if (sw == 3)
  {
    pop_size = popsize + new_row ; //轉存解的數量   
    for(i=1; i<=pop_size; i++)
    {
      for (j=1; j<=nobj; j++)
      {
        temp_pop[i].objv[j] = mixed_pop[i].objv[j] ; //cost value的轉存
      }
    }
  }
  for (i=1; i<=pop_size; i++)
    for (j=1; j<=pop_size; j++)
      _front[i].member[j] = 0 ;
     
  front_rank = 1 ;
  w = 0 ;

  for (i=1; i<=pop_size; i++)
  {
    temp_pop[i].nnp = 0 ;
    temp_pop[i].nsp = 0 ;
    temp_pop[i].front_rank = 0 ;
    temp_pop[i].crowd_dist = 0 ;
    _front[i].number = 0 ;
  }
  // 自己先跟所有人比 然後下一個跟所有人比...
  for (i=1; i<=pop_size; i++)
  {
    for (j=1; j<=pop_size; j++)
    {
      less=0 ;
      equal=0 ; 
      more=0 ;
      for (k=1; k<=nobj; k++) 
      {
        if (temp_pop[i].objv[k] < temp_pop[j].objv[k])
        {
          less ++ ;
        }
        else if (temp_pop[i].objv[k] == temp_pop[j].objv[k])
        {
          equal ++ ;
        }
        else
        {
          more ++ ;
        }
      }	

      if( less == 0 && equal != nobj) 
      {
        // 大於等於 或 大於大於
        temp_pop[i].nnp ++ ;  //計數器加1
        temp_pop[i].np[ temp_pop[i].nnp ] = j ;  
      }
      else if (more == 0 && equal != nobj) 
      {
        // 小於等於 或 小於小於
        temp_pop[i].nsp ++ ;  //計數器加1
        temp_pop[i].sp[ temp_pop[i].nsp ] = j ;  //??
      }	
  
    }

    if ( temp_pop[i].nnp == 0 ) 
    {
      //如果沒有"等號配上一個大於"的情況，那保證必定是小於的情況
      w ++ ; //計數器加1 
      temp_pop[i].front_rank = 1 ;  //  說明第i個解的排名是1 
      _front[1].member[w] = i ; // 相同front的member中對應到temp_pop的那個位置的解
    }
    _front[1].number = w ; //記下目前有幾個 front 1 的個數 
  }

	
///////////////////// 找出第二個Front以後，含有等號的相對應的rank  /////////////////////{ 
	
  while ( _front[front_rank].member[1] != 0 )       
  {
    temp1 = new int[2*popsize+1] ; // 整數指針
    w= 0 ;
    for (i=1; i<=_front[front_rank].number; i++)
    //front_rank值於此副程式的開頭定義成 1，而_front[front_rank].number表示記下目前有幾個 front 1 的個數
    {
      if ( temp_pop[ _front[front_rank].member[i] ].sp[1] != 0 )
      {
        for (j=1; j<=temp_pop[ _front[front_rank].member[i] ].nsp; j++) // temp_pop[] 裡的
        //"當第 i 個解為front 1 ，且在i的編號之前出現含有一個等號配上一個小於的情況 "的計數器數值拿來用
        {
          temp_pop[ temp_pop[ _front[front_rank].member[i] ].sp[j] ].nnp = temp_pop[ temp_pop[ _front[front_rank].member[i] ].sp[j] ].nnp - 1 ;     
          if ( temp_pop[ temp_pop[ _front[front_rank].member[i] ].sp[j] ].nnp == 0 ) // 若該解沒有被任何解支配, 則分配Front rank
          {
            temp_pop[ temp_pop[ _front[front_rank].member[i] ].sp[j] ].front_rank = front_rank + 1 ;
            w ++ ; // 計數每個Front成員的個數
            temp1[w] = temp_pop[ _front[front_rank].member[i] ].sp[j] ; 
          }
        }
      }
    }
    front_rank = front_rank + 1 ; // front_rank
    for (j=1; j<= w; j++)
      _front[front_rank].member[j] = temp1[j] ;

    _front[front_rank].number = w ; 
    delete [] temp1 ;
    temp1 = 0 ;
  } 
  front_size = front_rank-1 ; // 排除front1後，統計該回合產生多少個Front的轉存變數
  k=0 ;
  for (i=1; i<=pop_size; i++)
    front_sort[i] = 0 ; //歸零

  for (i=1; i<=front_size; i++)
  {
    for (j=1; j<=_front[i].number; j++)
    {
      k ++ ;
      front_sort[k] = _front[i].member[j] ;
    }
  }
  for (i=1; i<=pop_size; i++)
    temp_pop[ front_sort[i] ].sorted_front = i ;

  for (front_rank =1 ; front_rank<=front_size; front_rank++)
  {
    for (m=1; m<=nobj; m++)
    {
      temp_sort = new int[2*popsize+1] ;
      sort_on_obj(front_rank , m);  //按照 cost value 排序
      
      _fmin = temp_pop[ temp_sort[1] ].objv[m] ;  //cost value 的最小值為 front 1 
      _fmax = temp_pop[ temp_sort[ _front[front_rank].number ] ].objv[m] ; //cost value 的最大值為最後一名的 front  	
      
      temp_pop[ temp_sort[1] ].dist[m] = INF ;  // front 1 的 crowding distance 定為無窮大
      temp_pop[ temp_sort[ _front[front_rank].number ] ].dist[m] = INF ; // 最後一名的 front 的 crowding distance 定為無窮大

      for (j=2; j<=(_front[ front_rank].number-1); j++)
      {
        next_obj = temp_pop[ temp_sort[j+1] ].objv[m] ;
        pre_obj = temp_pop[ temp_sort[j-1] ].objv[m] ;
        
          if ( fabs(_fmax-_fmin) <= 1E-8 )
          temp_pop[ temp_sort[j] ].dist[m] = INF ;
          else
          temp_pop[ temp_sort[j] ].dist[m] = (next_obj-pre_obj)/(_fmax-_fmin) ;  //??
      }   
      ////////  釋放 temp1 的記憶體  ///////{
      delete [] temp_sort ;
      temp_sort = 0 ;
      ////////  釋放 temp1 的記憶體  ///////}
    }  
/////////  累加一組解的所有cost value 的 crowding distance ////////{   
    for (i=1; i<= _front[front_rank].number; i++)
    {
      _dist = 0;
      for (m=1; m<=nobj; m++)
        _dist = _dist + temp_pop[ _front[front_rank].member[i] ].dist[m] ; //累加一組解的所有cost value 的 crowding distance  
      temp_pop[ _front[front_rank].member[i] ].crowd_dist = _dist ; //轉存累加後的值
    }
//////////  累加一組解的所有cost value 的 crowding distance //////////}
  } 	

  if ( sw == 1)
  {
    for (i=1; i<=pop_size; i++)
    {
    parent_pop[i].nnp = temp_pop[i].nnp ; // 含有一個等號配上一個大於的計數器的轉存
    parent_pop[i].nsp = temp_pop[i].nsp ; // 含有一個等號配上一個小於的計數器的轉存
    parent_pop[i].sorted_front = temp_pop[i].sorted_front ; //將solutions按照Front rank排序前的轉存
    parent_pop[i].front_rank = temp_pop[i].front_rank ; //轉存 rank
    parent_pop[i].crowd_dist = temp_pop[i].crowd_dist ; //轉存 crowding distance ，用在sort_on_crowd的函式
    
      for (j=1; j<=nobj; j++)
      {
      parent_pop[i].dist[j] = temp_pop[i].dist[j]; //轉存累加的 crowding distance
      }
      
      for (j=1; j<= temp_pop[i].nsp; j++)
      {
      parent_pop[i].sp[j] = temp_pop[i].sp[j] ;  // 轉存含有一個等號配上一個小於的情況的數量解??
      }
    }
  }
  else if (sw == 3)
  {
    for (i=1; i<=pop_size; i++)
    {
    mixed_pop[i].nnp = temp_pop[i].nnp ;// 含有一個等號配上一個大於的計數器的轉存
    mixed_pop[i].nsp = temp_pop[i].nsp ; // 含有一個等號配上一個小於的計數器的轉存
    mixed_pop[i].sorted_front = temp_pop[i].sorted_front ;//將solutions按照Front rank排序前的轉存
    mixed_pop[i].front_rank = temp_pop[i].front_rank ;//轉存 rank
    mixed_pop[i].crowd_dist = temp_pop[i].crowd_dist ;//轉存 crowding distance ，用在sort_on_crowd的函式
    
      for (j=1; j<=nobj; j++)
      {
       mixed_pop[i].dist[j] = temp_pop[i].dist[j] ;//轉存累加的 crowding distance
      }
    
      for (j=1; j<= temp_pop[i].nsp; j++)
      {
       mixed_pop[i].sp[j] = temp_pop[i].sp[j] ;// 轉存含有一個等號配上一個小於的情況的數量解??
      }
    }
  }
//////////////////  子代的轉存  ////////////////////}   
}

void sort_on_crowd(int sw , int gg) 
{
  int pass,i,j,k,s,h=0;
  double hold;
  int pop_size=0 ;

  if (sw == 1)
  {
  pop_size = popsize ;

    for (i=1; i<=pop_size; i++)
    temp_pop[i].crowd_dist = parent_pop[i].crowd_dist ;
  }
  else if (sw == 3)
  {
  pop_size = popsize + new_row ;
  
  for (i=1; i<=pop_size; i++)
  temp_pop[i].crowd_dist = mixed_pop[i].crowd_dist ;
  }

  for (i=1; i<=pop_size; i++)
    sorted_crowding[i] = 0 ;  //暫存變數的歸零
  
  for (k=1; k<=front_size; k++)
  {
    for (j=1; j<=popsize + new_row; j++)  //暫存變數的歸零
    {
    temp_crowd[j][1] = 0. ;
    temp_crowd[j][2] = 0. ;
    }
    
    for (j=1; j<=_front[k].number; j++)
    {
    temp_crowd[j][1] = _front[k].member[j] ;  //轉存解的個數
    temp_crowd[j][2] = temp_pop[ _front[k].member[j] ].crowd_dist ; //轉存crodwing distance
    }
    ///////////////////////  排序，順序正確即保存  /////////////////////{ 
    for (pass=1; pass<= _front[k].number-1; pass++)
    {
      for (i=1; i<= _front[k].number-1; i++)
      {
        if (temp_crowd[i][2] < temp_crowd[i+1][2])
        {
          for (j=1; j<=2; j++)
          {
           //這裡使用變數的交換儲存技巧
	hold = temp_crowd[i][j] ;
	temp_crowd[i][j] = temp_crowd[i+1][j] ;
	temp_crowd[i+1][j] = hold ;
          }
        }
      }
    }
    ///////////////////////  排序，順序正確即保存  /////////////////////}

    for (s=1; s<=_front[k].number; s++)
    {
     h ++ ;
     sorted_crowding[h] = temp_crowd[s][1] ; //將排序好的結果轉存  (尚未對應到popsize)
    }	
  }
  
   for (i=1; i<=pop_size; i++)  
   temp_pop[ sorted_crowding[i] ].sorted_crowd =  i ; //設定的數量以符合接下來的轉存，因為親代的數量與子代不同。
   
//////////////////  親代的轉存  ////////////////////{   
   if (sw == 1)
   {
   for (i=1; i<=pop_size; i++)
   parent_pop[i].sorted_crowd = temp_pop[i].sorted_crowd ;
   }
//////////////////  親代的轉存  ////////////////////}   

//////////////////  子代的轉存  ////////////////////{   
   else if (sw == 3)
   {
     for (i=1; i<=pop_size; i++)
     {
     mixed_pop[i].sorted_crowd = temp_pop[i].sorted_crowd ;
     }
    sort_on_rank_crowd(pop_size) ; //做排列子代解的多樣性 
   }
//////////////////  子代的轉存  ////////////////////}
}

void sort_on_rank_crowd(int pop_size) 
{
// pop_size 是從sort_on_crowd內傳過來的
int pass,i,j;
double hold;
double temp_obj[2*popsize+1][nobj+2] ; //因為子代解至少會比親代多，所以創造2倍的空間

  for (i=1 ; i<= pop_size; i++) //注意這裡只存"前40個解"
  {
  temp_obj[i][1] = i ; //轉存解的個數
  temp_obj[i][2] = mixed_pop[i].sorted_crowd ;//轉存crodwing distance	
  }
 ///////////////////////  排序，順序正確即保存  /////////////////////{   
  for (pass=1 ; pass<= pop_size-1; pass++)
  {
    for( i=1; i<= pop_size-1; i++)
    {
      if (temp_obj[i][2] > temp_obj[i+1][2] )
      {
        for (j=1; j<=2; j++)
        {
         hold = temp_obj[i][j] ;
         temp_obj[i][j] = temp_obj[i+1][j] ;
         temp_obj[i+1][j] = hold;
        }
      }
    }
  }
 ///////////////////////  排序，順序正確即保存  /////////////////////}
 
  for (i=1 ; i<= pop_size; i++)
  {
  mixed_pop_sorted[i] = temp_obj[i][1] ; //排序完畢，轉存 
  }

}

void report_best(int c_gen, int sw )
{
  double sorted_temp[2*popsize+1][nobj+1] ;
  int i,pass ;  
  for (int j=1; j<=nobj; j++)
  {
    int w = 1 ;
    for (i=1; i<=popsize;i++)
    {  
      sorted_temp[i][j] = parent_pop[i].objv[j] ;  //轉存 cost value 
    } 
    for (pass=2; pass<=popsize-1; pass++)
    {
      if (sorted_temp[w][j] > sorted_temp[pass][j])   // 為何是比第一個costvalue??
      {
        best_solution[j] = pass ;
        w=pass;
      } 
    }
    if(w==1)
    {
      best_solution[j]=w; 
    }
  }
}


void update_pop(double _cg , double _mg , int kk )
{
  int _tour1 , _tour2 , jj  ;
  for (jj=1; jj<=popsize; jj=jj+2)
  {	
    do
    {
      _tour1 = tour_selection() ;  // 使用競爭法
      _tour2 = tour_selection() ;  // 使用競爭法
    }
    while (_tour1==_tour2);   
    crossover(jj , _tour1, _tour2) ; //交配
  }	

  mutation();  //突變		
}

int tour_selection() 
{
  int *tour_cand , *cand_rank;
  double *cand_dis ;
  int tour_member , j ;
  int tour_size = 2 ;
  tour_cand = new int[tour_size+1] ;
  cand_rank = new int[popsize+1] ;
  cand_dis = new double[popsize+1] ;
  for (j=1; j<=tour_size; j++)
  {
    tour_cand[j] = random(popsize) + 1 ;
    cand_rank[j] = parent_pop[ tour_cand[j] ].front_rank ;
    cand_dis[j] = parent_pop[ tour_cand[j] ].crowd_dist ;
  }
  
  if (cand_rank[1] == cand_rank[2])
  {
    if (cand_dis[1] >= cand_dis[2])
      tour_member = tour_cand[1] ;
    else
      tour_member = tour_cand[2] ;
  }
  else
  {
    if (cand_rank[1] < cand_rank[2])
      tour_member = tour_cand[1];
    else
      tour_member = tour_cand[2];
  }

  delete [] tour_cand;
  delete [] cand_rank;
  delete [] cand_dis ;
  tour_cand = 0 ;
  cand_rank = 0 ;
  cand_dis = 0 ;

  return tour_member;
}

void crossover (int kk , int _tour1 , int _tour2)
{
  int i;
  float rnd,par1,par2,chld1,chld2,betaq,beta,alpha;
  float y1,y2,y_max,y_min,expp;  
  rnd=randn(0,1);
  if(rnd<=_pcross)
  {
    for(i=1; i<=real_mem; i++)
    {
      par1=parent_pop[_tour1].xreal[i];
      par2=parent_pop[_tour2].xreal[i];
      y_max=max_variable(i);
      y_min=min_variable(i);
      
      rnd=randn(0,1);
      
      if(rnd<=0.5)
      {
        if(fabs(par1 - par2) > 0.000001) 
        {
          if(par2 > par1)
          {
            y2 = par2;
            y1 = par1;
          }
          else
          {
            y2 = par1;
            y1 = par2;
          }
                                           
          if((y1 - y_min) > (y_max - y2))     /*Find beta value*/
          {
            beta = 1 + (2*(y_max - y2)/(y2 - y1));
          }
          else
          {
            beta = 1 + (2*(y1-y_min)/(y2-y1));
          }                            
          expp = eta_c + 1.0; 
          beta = 1.0/beta; 
          alpha = 2.0 - pow(beta,expp);                              
          if (alpha < 0.0) 
          {
            alpha=0.000001;
          }                               
          rnd=randn(0,1);                                  
          if (rnd <= 1.0/alpha)
          {
            alpha = alpha*rnd;
            expp = 1.0/(eta_c+1.0);
            betaq = pow(alpha,expp);
          }
          else
          {
            alpha = alpha*rnd;
            alpha = 1.0/(2.0-alpha);
            expp = 1.0/(eta_c+1.0);
            if (alpha < 0.0) 
              alpha=0.000001;
            if (alpha < 0.0) 
            {
             exit(-1);
            }                                       
            betaq = pow(alpha,expp);                                          
          }                                           
          chld1 = 0.5*((y1+y2) - betaq*(y2-y1));
          chld2 = 0.5*((y1+y2) + betaq*(y2-y1));
       }
       else
       {                                  
        betaq = 1.0;
        y1 = par1; y2 = par2;
                                        
         /*Generation two children*/
         chld1 = 0.5*((y1+y2) - betaq*(y2-y1));
         chld2 =  0.5*((y1+y2) + betaq*(y2-y1));
       }
      if (chld1 < y_min) chld1 = y_min;
      if (chld1 > y_max) chld1 = y_max;
      if (chld2 < y_min) chld2 = y_min;
      if (chld2 > y_max) chld2 = y_max;                           
    }
    else
    {
      chld1 = par1;
      chld2 = par2;
    }
    child_pop[kk].xreal[i]=chld1;
    child_pop[kk+1].xreal[i]=chld2;
  }		
 }
}

void mutation()
{
  int i,j;
  float rnd,delta,indi,deltaq;
  float y,y_min,y_max,val,xy;
  for (i=1; i<=popsize; i++)
  {
    for(j=1; j<=real_mem; j++)
    {
    rnd=randn(0,1);
     if(rnd<= _pmutat) // for each variable find whether to do mutation or not
     {
      y=child_pop[i].xreal[j];
      y_max=max_variable(j);
      y_min=min_variable(j);                
       if (y>y_min)
       {
         if((y-y_min) < (y_max-y))
          delta = (y - y_min)/(y_max - y_min);
                                
         else
         delta = (y_max - y)/(y_max - y_min);
                                
         rnd=randn(0,1);
         indi=1.0/(eta_m+1.0);
                                
         if (rnd<=0.5)
         {
          xy = 1.0-delta;
          val = 2*rnd+(1-2*rnd)*(pow(xy,(eta_m+1)));
          deltaq =  pow(val,indi) - 1.0;
         }
         else
         {
          xy = 1.0-delta;
          val = 2.0*(1.0-rnd)+2.0*(rnd-0.5)*(pow(xy,(eta_m+1)));
          deltaq = 1.0 - (pow(val,indi));
         }
          
          y = y + deltaq * (y_max-y_min);
                                 
          if (y < y_min) y=y_min;
          if (y > y_max) y=y_max;
                                 
          child_pop[i].xreal[j]=y;
        }
        else
        {
         xy=randn(0,1);
         child_pop[i].xreal[j]=xy*(y_max-y_min)+y_min;
        }
      }
    }
  }
}


void child_pop_cluster(int i ,int step,int c_gen)
{
    int ii,iii,jj,jjj;
    for(jj=1; jj<=_in_clu; jj++)
	  {
	    for(jjj=1; jjj<=_in_varl; jjj++)
		  {
		    fuzzy_real[i][(jj-1)*_mem_length + jjj*2 - 1] = child_pop[i].xreal[(jj-1)*_mem_length + jjj*2 - 1] ; // 中心點轉存
        fuzzy_real[i][(jj-1)*_mem_length + jjj*2] = child_pop[i].xreal[(jj-1)*_mem_length + jjj*2] ; // 寬度轉存
        if(fuzzy_real[i][(jj-1)*_mem_length+jjj*2-1] < min_m[i][jjj])
        {
          min_m[i][jjj]=fuzzy_real[i][(jj-1)*_mem_length+jjj*2-1];  //如果比最大的最小值(min_m = 1)還小，就儲存最小中心點
        }
        
        if(fuzzy_real[i][(jj-1)*_mem_length+jjj*2-1] > max_m[i][jjj])
        {
          max_m[i][jjj]=fuzzy_real[i][(jj-1)*_mem_length+jjj*2-1];  //如果比最小的最大值(max_m = 1)還大，就儲存最大中心點，
        }   
		  }
		  
		  // fuzzy_real[i][jj*_mem_length-_out_varl+1] =  child_pop[i].xreal[jj*_mem_length-_out_varl+1] * left_wheel_speed ; // 左輪
      // fuzzy_real[i][jj*_mem_length-_out_varl+2] =  child_pop[i].xreal[jj*_mem_length-_out_varl+2] * right_wheel_speed ;  // 右輪				
      for(int k=1; k<=9; k++) // 後件部轉存
      {
        fuzzy_real[i][jj*_mem_length-_conq_varl + k] =  child_pop[i].xreal[jj*_mem_length - _conq_varl + k]; 
      }
	  }
}

void merge() 
{
int i, w;	
for (i=1; i<=popsize; i++)
  mixed_pop[i] = parent_pop[i] ; 

for (i=1, w=(popsize+1); i<=new_row; i++, w++)
  mixed_pop[w] = child_pop[i] ;
}

void update_newparent()
{
	int i,j,jj ;
	for(i=1;i<=popsize;i++)
		parent_pop[i] = mixed_pop[ mixed_pop_sorted[i] ] ;
	//child_pop initial
	for(i=1; i<=popsize + new_row; i++)
	{
		child_pop[i].front_rank = 0 ;
		child_pop[i].crowd_dist = 0 ;
				
		for (j=1; j<=_max_rule*_mem_length; j++)
			child_pop[i].xreal[j] = 0 ;

		for (j=1; j<=nobj; j++)
			child_pop[i].objv[j] = 0 ;
	}

	//temp_child initial
	for(i=1; i<=popsize+new_row; i++)
	{
		temp_pop[i].front_rank = 0 ;
		temp_pop[i].crowd_dist = 0 ;
		
		for (j=1; j<=_max_rule*_mem_length; j++)
			temp_pop[i].xreal[j] = 0 ;

		for (j=1; j<=nobj; j++)
			temp_pop[i].objv[j] = 0 ;
	}
	
	//mixed_pop initial
	for(i=1; i<=popsize+new_row; i++)
	{
		mixed_pop[i].front_rank = 0 ;
		mixed_pop[i].crowd_dist = 0 ;
		
		for (j=1; j<=_max_rule*_mem_length; j++)
			mixed_pop[i].xreal[j] = 0 ;

		for (j=1; j<=nobj; j++)
			mixed_pop[i].objv[j] = 0 ;
	}	
}

void save_final_pareto_front()
{
  int i,j ;
  FILE *fpt_init; 
  fpt_init = fopen("..\\..\\generated_file\\final_pareto_front00.txt","a") ;
  for (i=1; i<=popsize; i++) 
  {  
    for(j=1; j<=nobj;j++)  
    {
      fprintf(fpt_init,"%f\t", parent_pop[i].objv[j]) ;  //cost value
    }
   fprintf(fpt_init,"%d\r\n", parent_pop[i].front_rank); //排名
  } 
  fclose(fpt_init) ;
  FILE *fpt_init01;
  fpt_init01 = fopen("..\\..\\generated_file\\final_pareto_front.txt","a") ; 
  for (i=1; i<=popsize; i++)
  {
    for(j=1; j<=nobj;j++) 
    {
       fprintf(fpt_init01,"%f\t", parent_pop[i].objv[j]) ;  //cost value		
    }
    fprintf(fpt_init01, "\r\n") ;
  }
  fclose(fpt_init01) ;

}

void save_final_obj(int j) 
{	
  if (j==1)
  {	
    int jj , jjj ; 
    FILE *f_obj ;
    if( (f_obj=fopen("..\\..\\generated_file\\w.txt","a"))==NULL) exit(1) ;
    {
      for(jj=1; jj<=_rule_number; jj++)
      { 
        for(jjj=1; jjj<=_in_varl; jjj++)
        {
          fprintf(f_obj,"%.10f\t", parent_pop[best_solution[j]].xreal[(jj-1)*_mem_length+jjj*2-1] ); // 中心點
          fprintf(f_obj,"%.10f\t", parent_pop[best_solution[j]].xreal[(jj-1)*_mem_length+jjj*2] ) ; // 寬度
        }
        // fprintf(f_obj,"%.10f   \t",  parent_pop[best_solution[j]].xreal[jj*_mem_length-_out_varl+1] * left_wheel_speed ) ; // 左輪
        // fprintf(f_obj,"%.10f   \t",  parent_pop[best_solution[j]].xreal[jj*_mem_length-_out_varl+2] * right_wheel_speed  ) ;  // 右輪  
        for(int k = 1; k <= 9; k++)
        {
          fprintf(f_obj,"%.10f   \t",  parent_pop[best_solution[j]].xreal[jj*_mem_length - _conq_varl + k]) ;  
        }

      }
    }
    fprintf(f_obj,"\n") ;
    fclose(f_obj) ;
  }
  if (j==1)
  {	
    int jj ; 
    FILE *f_obj ;
    if( (f_obj=fopen("..\\..\\generated_file\\cost_w.txt","a"))==NULL) exit(1) ;
    {
      for(jj=1; jj<=nobj;jj++)  
      {
          fprintf(f_obj,"%f\t", parent_pop[best_solution[j]].objv[jj]) ;  //cost value
      }
    }
    fprintf(f_obj,"\n") ;
    fclose(f_obj) ;
  }
  if(j==2)
  {
    int jj , jjj ; 
    FILE *f_obj ;
    if( (f_obj=fopen("..\\..\\generated_file\\s.txt","a"))==NULL) exit(1) ;
    {
      for(jj=1; jj<=_rule_number; jj++)
      { 
        for(jjj=1; jjj<=_in_varl; jjj++)
        {
          fprintf(f_obj,"%.10f\t", parent_pop[best_solution[j]].xreal[(jj-1)*_mem_length+jjj*2-1] ); // 中心點
          fprintf(f_obj,"%.10f\t", parent_pop[best_solution[j]].xreal[(jj-1)*_mem_length+jjj*2] ) ; // 寬度
        }
        // fprintf(f_obj,"%.10f   \t",  parent_pop[best_solution[j]].xreal[jj*_mem_length-_out_varl+1] * left_wheel_speed ) ; // 左輪
        // fprintf(f_obj,"%.10f   \t",  parent_pop[best_solution[j]].xreal[jj*_mem_length-_out_varl+2] * right_wheel_speed  ) ;  // 右輪
        for(int k = 1; k <= 9; k++)
        {
          fprintf(f_obj,"%.10f   \t",  parent_pop[best_solution[j]].xreal[jj*_mem_length - _conq_varl + k]) ;  
        }  
      }
    }
    fprintf(f_obj,"\n") ;
    fclose(f_obj) ;
  }
  if (j==2)
  {	
    int jj ; 
    FILE *f_obj ;
    if( (f_obj=fopen("..\\..\\generated_file\\cost_s.txt","a"))==NULL) exit(1) ;
    {
      for(jj=1; jj<=nobj;jj++)  
      {
          fprintf(f_obj,"%f\t", parent_pop[best_solution[j]].objv[jj]) ;  //cost value
      }
    }
    fprintf(f_obj,"\n") ;
    fclose(f_obj) ;
  }

  if(j==3)
  {
    int jj , jjj ; 
    FILE *f_obj ;
    if( (f_obj=fopen("..\\..\\generated_file\\i.txt","a"))==NULL) exit(1) ;
    {
      for(jj=1; jj<=_rule_number; jj++)
      { 
        for(jjj=1; jjj<=_in_varl; jjj++)
        {
        fprintf(f_obj,"%.10f\t", parent_pop[best_solution[j]].xreal[(jj-1)*_mem_length+jjj*2-1] ); // 中心點
        fprintf(f_obj,"%.10f\t", parent_pop[best_solution[j]].xreal[(jj-1)*_mem_length+jjj*2] ) ; // 寬度
        }
        // fprintf(f_obj,"%.10f   \t",  parent_pop[best_solution[j]].xreal[jj*_mem_length-_out_varl+1] * left_wheel_speed ) ; // 左輪
        // fprintf(f_obj,"%.10f   \t",  parent_pop[best_solution[j]].xreal[jj*_mem_length-_out_varl+2] * right_wheel_speed  ) ;  // 右輪  
        for(int k = 1; k <= 9; k++)
        {
          fprintf(f_obj,"%.10f   \t",  parent_pop[best_solution[j]].xreal[jj*_mem_length - _conq_varl + k]) ;  
        }
      }
    }
    fprintf(f_obj,"\n") ;
    fclose(f_obj) ;
  }

  if (j==3)
  {	
  int jj ; 
  FILE *f_obj ;

    if( (f_obj=fopen("..\\..\\generated_file\\cost_i.txt","a"))==NULL) exit(1) ;
    {
      for(jj=1; jj<=nobj;jj++)  
      {
          fprintf(f_obj,"%f\t", parent_pop[best_solution[j]].objv[jj]) ;  //cost value
      }
    }
  fprintf(f_obj,"\n") ;
  fclose(f_obj) ;
    }

}

int end_iteration_and_save_gen(int load_gen)
{
  load_gen = load_gen+1;
  FILE *fnoise1;  
  fnoise1 = fopen("gen_record.txt","w");
  fprintf(fnoise1,"%d\n", load_gen) ;  
  fclose(fnoise1);
  printf("end_iteration_and_save_gen=%d\n", load_gen);
  fflush(stdout);
  return load_gen;
}

void load_save_parent_rule_and_relative_data()
{
  FILE *fnoise1;  
  fnoise1 = fopen("parent_rank.txt","r"); 	
  for (int i=1; i<=popsize; i++)
  {
    fscanf(fnoise1,"%d\t",&parent_rank_load[i]);	        
    fscanf(fnoise1, "\n") ;
  }
  fclose(fnoise1);
  for(int i=1; i<=popsize; i++ )
  {   
    parent_pop[i].front_rank=parent_rank_load[i];          
  }
  FILE *fnoise2;
  fnoise2 = fopen("parent_obj.txt","r"); 	 
  for (int i=1; i<=popsize; i++) 
  {   
    for(int jj=1; jj<=nobj; jj++)
    {     
      fscanf(fnoise2,"%Lf\t",&parent_obj_load[i][jj]);	        			
    }    
    fscanf(fnoise2, "\n") ;
  }   
  fclose(fnoise2);
  for(int i=1; i<=popsize; i++ )
  {
    for(int jj=1; jj<=nobj; jj++)
    {
      printf("Reading objv[%d][%d]: %Lf\n", i, jj, parent_obj_load[i][jj]); // wei Debug
      parent_pop[i].objv[jj]=parent_obj_load[i][jj];        
    }   
  }
  FILE *fnoise3;
  fnoise3 = fopen("parent_rule.txt","r"); 	
  for (int i=1; i<=popsize; i++) 
  {   
    for(int jj=1; jj<=real_mem; jj++)
    {     
        fscanf(fnoise3,"%Lf\t",&parent_rule_load[i][jj]);	        			
    }    
    fscanf(fnoise3, "\n") ;
  }   
  fclose(fnoise3);

  for(int i=1; i<=popsize; i++ )
  {
    for(int jj=1; jj<=_rule_number; jj++)
    {
        for(int jjj=1; jjj<=_in_varl; jjj++)
        {        
          parent_pop[i].xreal[(jj-1)*_mem_length+jjj*2-1] = parent_rule_load[i][(jj-1)*_mem_length+jjj*2-1] ; // 中心點
          parent_pop[i].xreal[(jj-1)*_mem_length+jjj*2] = parent_rule_load[i][(jj-1)*_mem_length+jjj*2] ; // 寬度
              
        }
        // parent_pop[i].xreal[jj*_mem_length-_out_varl+1] = parent_rule_load[i][jj*_mem_length-_out_varl+1]; // 左輪
        // parent_pop[i].xreal[jj*_mem_length-_out_varl+2] = parent_rule_load[i][jj*_mem_length-_out_varl+2];  // 右輪
        for(int k=1; k<=9; k++) // 後件部參數讀檔
        {
          parent_pop[i].xreal[jj*_mem_length - _conq_varl + k] = parent_rule_load[i][jj*_mem_length - _conq_varl + k];
        }

    }   
  }
  fflush(stdout);
}

double randn(double nmin, double nmax)
{
	double rinv;
	rinv=nmax-nmin;
	return ( (double)(rand()%(RAND_MAX+1)) / (double)RAND_MAX )*rinv+nmin;
}

void gauss_data_limit(int kk)
{
int j,jj;
  for (jj=1; jj<=_rule_number; jj++)
  {
    for(j=1;j<=_in_varl;j++)
    {
	   if (child_pop[kk].xreal[(jj-1)*_mem_length+j*2-1] < _center_min)
      {
       child_pop[kk].xreal[(jj-1)*_mem_length+j*2-1] = _center_min ; // 限制中心點最小值
      }
 
      if (child_pop[kk].xreal[(jj-1)*_mem_length+j*2-1] > _center_max)
      {
         child_pop[kk].xreal[(jj-1)*_mem_length+j*2-1] = _center_max ; // 限制中心點最大值   
      }
	  
		
      if (child_pop[kk].xreal[(jj-1)*_mem_length+j*2] < _range_min)
      {
       child_pop[kk].xreal[(jj-1)*_mem_length+j*2] = _range_min ; // 限制寬度最小值
      }
 
      if (child_pop[kk].xreal[(jj-1)*_mem_length+j*2] > _range_max)
      {
         child_pop[kk].xreal[(jj-1)*_mem_length+j*2] = _range_max ; // 限制寬度最大值   
      }
    }		
  }
}

void save_temp_obj(int j) 
{
	
if (j==1)
{	
int jj , jjj ; 
FILE *f_obj ;

  if( (f_obj=fopen("..\\..\\generated_file\\temp_check\\temp_w.txt","w"))==NULL) exit(1) ;
 // if( (f_obj=fopen("fuzzy_real_temp_obj01_best.txt","w"))==NULL) exit(1) ;
  {
    for(jj=1; jj<=_rule_number; jj++)
    { 
      for(jjj=1; jjj<=_in_varl; jjj++)
      {
      fprintf(f_obj,"%.10f\t", parent_pop[best_solution[j]].xreal[(jj-1)*_mem_length+jjj*2-1] ); // 中心點
      fprintf(f_obj,"%.10f\t", parent_pop[best_solution[j]].xreal[(jj-1)*_mem_length+jjj*2] ) ; // 寬度
      }
    // fprintf(f_obj,"%.10f   \t",  parent_pop[best_solution[j]].xreal[jj*_mem_length-_out_varl+1] * left_wheel_speed ) ; // 左輪
    // fprintf(f_obj,"%.10f   \t",  parent_pop[best_solution[j]].xreal[jj*_mem_length-_out_varl+2] * right_wheel_speed  ) ;  // 右輪  
    for(int k = 1; k <= 9; k++)
    {
      fprintf(f_obj,"%.10f   \t",  parent_pop[best_solution[j]].xreal[jj*_mem_length - _conq_varl + k]) ;  
    }
    }
  }
fprintf(f_obj,"\n") ;
fclose(f_obj) ;
}

if (j==1)
{	
int jj ; 
FILE *f_obj ;

  if( (f_obj=fopen("..\\..\\generated_file\\temp_check\\cost_temp_w.txt","w"))==NULL) exit(1) ;
  {
     for(jj=1; jj<=nobj;jj++)  
     {
         fprintf(f_obj,"%f\t", parent_pop[best_solution[j]].objv[jj]) ;  //cost value
     }
  }
fprintf(f_obj,"\n") ;
fclose(f_obj) ;
}

if(j==2)
{
int jj , jjj ; 

FILE *f_obj ;

  if( (f_obj=fopen("..\\..\\generated_file\\temp_check\\temp_s.txt","w"))==NULL) exit(1) ;
  {
    for(jj=1; jj<=_rule_number; jj++)
    { 
      for(jjj=1; jjj<=_in_varl; jjj++)
      {
      fprintf(f_obj,"%.10f\t", parent_pop[best_solution[j]].xreal[(jj-1)*_mem_length+jjj*2-1] ); // 中心點
      fprintf(f_obj,"%.10f\t", parent_pop[best_solution[j]].xreal[(jj-1)*_mem_length+jjj*2] ) ; // 寬度
      }
    // fprintf(f_obj,"%.10f   \t",  parent_pop[best_solution[j]].xreal[jj*_mem_length-_out_varl+1] * left_wheel_speed ) ; // 左輪
    // fprintf(f_obj,"%.10f   \t",  parent_pop[best_solution[j]].xreal[jj*_mem_length-_out_varl+2] * right_wheel_speed  ) ;  // 右輪  
    for(int k = 1; k <= 9; k++)
    {
      fprintf(f_obj,"%.10f   \t",  parent_pop[best_solution[j]].xreal[jj*_mem_length - _conq_varl + k]) ;  
    }
    }
  }
fprintf(f_obj,"\n") ;
fclose(f_obj) ;
}

if (j==2)
{	
int jj ; 
FILE *f_obj ;

  if( (f_obj=fopen("..\\..\\generated_file\\temp_check\\cost_temp_s.txt","w"))==NULL) exit(1) ;
  {
     for(jj=1; jj<=nobj;jj++)  
     {
         fprintf(f_obj,"%f\t", parent_pop[best_solution[j]].objv[jj]) ;  //cost value
     }
  }
fprintf(f_obj,"\n") ;
fclose(f_obj) ;
}

if(j==3)
{
int jj , jjj ; 
FILE *f_obj ;
 
  if( (f_obj=fopen("..\\..\\generated_file\\temp_check\\temp_i.txt","w"))==NULL) exit(1) ;
  {
    for(jj=1; jj<=_rule_number; jj++)
    { 
      for(jjj=1; jjj<=_in_varl; jjj++)
      {
      fprintf(f_obj,"%.10f\t", parent_pop[best_solution[j]].xreal[(jj-1)*_mem_length+jjj*2-1] ); // 中心點
      fprintf(f_obj,"%.10f\t", parent_pop[best_solution[j]].xreal[(jj-1)*_mem_length+jjj*2] ) ; // 寬度
      }
    // fprintf(f_obj,"%.10f   \t",  parent_pop[best_solution[j]].xreal[jj*_mem_length-_out_varl+1] * left_wheel_speed ) ; // 左輪
    // fprintf(f_obj,"%.10f   \t",  parent_pop[best_solution[j]].xreal[jj*_mem_length-_out_varl+2] * right_wheel_speed  ) ;  // 右輪
    for(int k = 1; k <= 9; k++)
    {
      fprintf(f_obj,"%.10f   \t",  parent_pop[best_solution[j]].xreal[jj*_mem_length - _conq_varl + k]) ;  
    }  
    }
  }
fprintf(f_obj,"\n") ;
fclose(f_obj) ;
}

if (j==3)
{	
int jj ; 
FILE *f_obj ;

  if( (f_obj=fopen("..\\..\\generated_file\\temp_check\\cost_temp_i.txt","w"))==NULL) exit(1) ;
  {
     for(jj=1; jj<=nobj;jj++)  
     {
         fprintf(f_obj,"%f\t", parent_pop[best_solution[j]].objv[jj]) ;  //cost value
     }
  }
fprintf(f_obj,"\n") ;
fclose(f_obj) ;
}
 
}

void report_pop ( FILE *fpt)
{
int i, j , jj , jjj ;

for (i=1; i<=popsize; i++)
{
  for (j=1; j<=nobj; j++)
  {
   fprintf(fpt,"%.10f\t", parent_pop[i].objv[j]);
  }
   fprintf(fpt," , \t");
		
  for(jj=1; jj<=_rule_number; jj++)
  {
    for(jjj=1; jjj<=_in_varl; jjj++)
    {
    fprintf(fpt,"%.10f\t", parent_pop[i].xreal[(jj-1)*_mem_length+jjj*2-1] ); // 中心點
    fprintf(fpt,"%.10f\t", parent_pop[i].xreal[(jj-1)*_mem_length+jjj*2] ) ; // 寬度
    }
    // fprintf(fpt,"%.10f   \t",  parent_pop[i].xreal[jj*_mem_length-_out_varl+1] * left_wheel_speed ) ; // 左輪
    // fprintf(fpt,"%.10f   \t",  parent_pop[i].xreal[jj*_mem_length-_out_varl+2] * right_wheel_speed  ) ;  // 右輪  
    for(int k=1; k<=9; k++)
    {
      fprintf(fpt,"%.10f\t",  parent_pop[i].xreal[jj*_mem_length - _conq_varl + k]) ; // wei debug
    }	
     	

   // fuzzy_real[i][jj*_mem_length-_out_varl+1] = parent_pop[i].xreal[jj*_mem_length-_out_varl+1] * left_wheel_speed ; // 左輪
  // fuzzy_real[i][jj*_mem_length-_out_varl+2] =  parent_pop[i].xreal[jj*_mem_length-_out_varl+2] * right_wheel_speed ;  // 右輪
  // 後件部轉存
  // for(int k=1; k<=9; k++)
  // {
  //   fuzzy_real[i][jj*_mem_length-_conq_varl + k] =  parent_pop[i].xreal[jj*_mem_length - _conq_varl + k]; 
  // }	
  }

fprintf(fpt," , \t");
fprintf(fpt,"%d\t", parent_pop[i].front_rank);
fprintf(fpt,"%.10f\n", parent_pop[i].crowd_dist);
}
return;
}

void sort_on_obj(int front_rank , int m) 
{
//這裡的 m 是 nobj 。

int pass,i,j;
double hold;
double temp_obj[2*popsize+1][nobj+2] ;

  for (i=1 ; i<= _front[front_rank].number; i++) 
  {
   temp_obj[i][1] = _front[front_rank].member[i] ; // 轉存第幾個解
   temp_obj[i][2] = temp_pop[ _front[front_rank].member[i] ].objv[m] ; //轉存cost value		
  }
  
 
///////////////////////  排序，順序正確即保存  /////////////////////{      
  for (pass=1 ; pass<= _front[ front_rank].number-1; pass++)
  {
    for( i=1; i<= _front[ front_rank].number-1; i++)
    {
      if (temp_obj[i][2] > temp_obj[i+1][2] )
      {
        for (j=1; j<=2; j++)
        {
        //這裡使用變數的交換儲存技巧
        hold = temp_obj[i][j] ;
        temp_obj[i][j] = temp_obj[i+1][j] ;
        temp_obj[i+1][j] = hold;
        }
      }
    }
  }
///////////////////////  排序，順序正確即保存  /////////////////////}

///////////////////////  將排好的順序重新轉存  ///////////////////{
  for (i=1 ; i<= _front[front_rank].number; i++)
  {
  temp_sort[i] = temp_obj[i][1] ;
  }
///////////////////////  將排好的順序重新轉存  ///////////////////}
}

void load_gen_data()
{
 
FILE *fnoise1;  

  fnoise1 = fopen("gen_record.txt","r");
   
    fscanf(fnoise1,"%d\n",&load_gen);	
  
  fclose(fnoise1);


printf("load_in=%d\n",load_gen);


fflush(stdout);
}