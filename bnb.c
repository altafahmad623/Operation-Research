#include<stdio.h>
#include<malloc.h>
#include<math.h>
#include <stdbool.h>
#include <stdlib.h>

int dual_simplex_method(double **T,int no_of_eqns,int no_of_variables,int *column,int *row){
    int flag=0,mp=0,i;
    double **A;
    A=(double**)malloc((no_of_eqns+1)*sizeof(double*));
    for(i=0;i<no_of_eqns+1;i++){
        A[i]=(double*)malloc((no_of_variables+1)*sizeof(double));
    }
    double min=100;
    for(i=0;i<no_of_eqns;i++){
      if(T[i][no_of_variables]<min&&T[i][no_of_variables]<0){
        flag=1;
        min=T[i][no_of_variables];
        mp=i;
      }
    }
    int most_negative_row=mp,pivot_column;
    double most_negative_value=min,pivot_value;
    if(flag==0){
        return 1;
    }else{
        int min_ratio_column=0;
        double min_ratio=10000;
        flag=0;
        for(i=0;i<no_of_variables;i++){
            if(fabs(T[no_of_eqns][i]/T[most_negative_row][i])<min_ratio&&T[most_negative_row][i]<0){
                min_ratio=fabs(T[no_of_eqns][i]/T[most_negative_row][i]);
                min_ratio_column=i;
            }
        }
        if(min_ratio==10000){
            printf("The solution is unbounded\n");
            return 2;
        }else{
            pivot_column=min_ratio_column;
            pivot_value=T[most_negative_row][min_ratio_column];
            int temp1=row[most_negative_row];
            row[most_negative_row]=column[min_ratio_column];
            column[min_ratio_column]=temp1;
            for(i=0;i<no_of_eqns+1;i++){
                int j;
                for(j=0;j<no_of_variables+1;j++){
                    if(i==most_negative_row&&j==min_ratio_column){
                        A[i][j]=1/T[i][j];
                    }
                    else if(i==most_negative_row&&j!=min_ratio_column){
                        A[i][j]=T[i][j]/pivot_value;
                    }
                    else if(j==min_ratio_column&&i!=most_negative_row){
                        A[i][j]=-T[i][j]/pivot_value;
                    }
                    else if(j!=min_ratio_column&&i!=most_negative_row){
                        double temp=T[i][j]*pivot_value-T[most_negative_row][j]*T[i][min_ratio_column];
                        A[i][j]=temp/((double)pivot_value);
                    }
                }
            }
        }
    }
    for(i=0;i<no_of_eqns+1;i++){
        int j;
        for(j=0;j<no_of_variables+1;j++){
                T[i][j]=A[i][j];
        }
    }
    printf("\nTable: \n");
    int k,l;
    for(l=0;l<no_of_variables;l++){
        if(column[l]/10==1){
            printf("X%d\t\t",(column[l]%10)+1);
        }else if(column[l]/10==2){
            printf("Z%d\t\t",(column[l]%10)+1);
        }
    }
    printf("1\n");
    for(k=0;k<no_of_eqns+1;k++){
      for(l=0;l<no_of_variables+1;l++){
        printf("%lf\t",T[k][l]);
      }
      if(k==no_of_eqns){
        printf("Z\n");
      }else{
          switch(row[k]/10){
            case 1:printf("X%d\n",(row[k]%10)+1);
                break;
            case 2:printf("Z%d\n",(row[k]%10)+1);
                break;
          }
      }
    }
    return 0;
}

void print_result(double **T,int no_of_eqns,int no_of_variables,int *column,int *row,int decision_variable_count){
    int i;
    printf("The maximum solution occurs at Z=%lf with\n",T[no_of_eqns][no_of_variables]);
    for(i=0;i<no_of_eqns;i++){
        if(row[i]/10==1&&row[i]%10<decision_variable_count){
            printf("X%d=%lf\t",(row[i]%10)+1,T[i][no_of_variables]);
        }
    }
    for(i=0;i<no_of_variables;i++){
        if(column[i]/10==1&&column[i]%10<decision_variable_count){
            printf("X%d=0\t",(column[i]%10)+1);
        }
    }
}


void print_result_phase_I(double **T,int no_of_eqns,int no_of_variables,int *column,int *row){
  int i;
  printf("The maximum solution occurs at Z=%lf with\n",T[no_of_eqns][no_of_variables]);
  for(i=0;i<no_of_eqns;i++){
    if(row[i]/10==2){
      printf("Z%d=%lf\t",(row[i]%10)+1,T[i][no_of_variables]);
    }
  }
  for(i=0;i<no_of_eqns;i++){
    if(column[i]/10==2){
      printf("Z%d=%lf\t",(column[i]%10)+1,0.00);
    }
  }
}

void print_result_phase_II(double **T,int no_of_eqns,int no_of_variables,int *column,int *row){
    int i;
    printf("The maximum solution occurs at Z=%lf with\n",T[no_of_eqns][no_of_variables]);
    for(i=0;i<no_of_eqns;i++){
        if(row[i]/10==1){
            printf("X%d=%lf\t",(row[i]%10)+1,T[i][no_of_variables]);
        }
    }
    for(i=0;i<no_of_variables;i++){
        if(column[i]/10==1){
            printf("X%d=%lf\t",(column[i]%10)+1,0.00);
        }
    }
}



int simplex_method(double **T,int no_of_eqns,int no_of_variables,int *column,int *row,int phase,int eqn_type[]){
  int flag=0,mp=0,i;
  double **A;
  int special1=0;
  A=(double**)malloc((no_of_eqns+1)*sizeof(double*));
  for(i=0;i<no_of_eqns+1;i++){
    A[i]=(double*)malloc((no_of_variables+1)*sizeof(double));
  }
  double min=100;
  for(i=0;i<no_of_variables;i++){
    if(T[no_of_eqns][i]<min&&T[no_of_eqns][i]<=0){
      flag=1;
      min=T[no_of_eqns][i];
      mp=i;
    }
  }
  int most_negative_column=mp,pivot_row;
  double most_negative_value=min,pivot_value=10000;
    repeat:
  if(flag==0){
    if(phase==1)
        print_result_phase_I(T,no_of_eqns,no_of_variables,column,row);
    else if(phase==2)
        print_result_phase_II(T,no_of_eqns,no_of_variables,column,row);
    return 1;
  }else{
    int min_ratio_row=0;
    double min_ratio=10000;
    int flag3=0;
    for(i=0;i<no_of_eqns;i++){
      if(T[i][no_of_variables]/T[i][most_negative_column]<min_ratio&&T[i][most_negative_column]>0){
        flag3=1;
        min_ratio=T[i][no_of_variables]/T[i][most_negative_column];
        min_ratio_row=i;
      }
    }


    if(phase==2){
        int j;
        for(j=0;j<no_of_eqns;j++){
            if(row[j]/10==2&&eqn_type[row[j]%10]!=1&&T[j][most_negative_column]<0){
                flag3=1;
                min_ratio=0;
                min_ratio_row=j;
                special1=1;
            }
        }
    }




    if(flag3==0){
        int k,temp=most_negative_column+1;
        for(k=temp;k<no_of_variables;k++){
            if(T[no_of_eqns][k]==0){
                most_negative_column=k;most_negative_value=0;
                goto repeat;
            }
        }
      printf("The solution is unbounded\n");
      return 2;
    }else{
      pivot_row=min_ratio_row;
      pivot_value=T[min_ratio_row][most_negative_column];
      int temp1=row[min_ratio_row];
      row[min_ratio_row]=column[most_negative_column];
      column[most_negative_column]=temp1;
      for(i=0;i<no_of_eqns+1;i++){
        int j;
        for(j=0;j<no_of_variables+1;j++){
            if(i==pivot_row&&j==most_negative_column){
                A[i][j]=1/T[i][j];
            }
            else if(i==pivot_row&&j!=most_negative_column){
                A[i][j]=T[i][j]/pivot_value;
            }else if(j==most_negative_column&&i!=pivot_row){
                A[i][j]=-T[i][j]/pivot_value;
            }else if(j!=most_negative_column&&i!=pivot_row){
                double temp=T[i][j]*pivot_value-T[pivot_row][j]*T[i][most_negative_column];
                A[i][j]=temp/((double)pivot_value);
            }
        }
      }
    }
  }
  for(i=0;i<no_of_eqns+1;i++){
    int j;
    for(j=0;j<no_of_variables+1;j++){
      T[i][j]=A[i][j];
    }
  }


  if(special1==1){
    int cc;
    for(i=0;i<no_of_variables;i++){
        if(column[i]/10==2&&eqn_type[column[i]%10]!=1){
            cc=i;break;
        }
    }
    for(i=0;i<=no_of_eqns;i++){
        T[cc][i]=0;
    }
  }


  printf("\nTable: \n");
  int k,l;
  for(l=0;l<no_of_variables;l++){
    if(column[l]/10==1){
      printf("X%d\t\t",(column[l]%10)+1);
    }else if(column[l]/10==2){
      printf("Z%d\t\t",(column[l]%10)+1);
    }
  }
  printf("1\n");
  for(k=0;k<no_of_eqns+1;k++){
    for(l=0;l<no_of_variables+1;l++){
      printf("%lf\t",T[k][l]);
    }
    if(k==no_of_eqns){
      printf("Z\n");
    }else{
      switch(row[k]/10){
      case 1:printf("X%d\n",(row[k]%10)+1);
    break;
      case 2:printf("Z%d\n",(row[k]%10)+1);
    break;
      }
    }
  }
  return 0;
}



int main(){
  int no_of_eqns=1,no_of_variables=1;
  printf("Enter the number of equations\n");
  scanf("%d",&no_of_eqns);
  printf("Enter the number of variables\n");
  scanf("%d",&no_of_variables);
  printf("\n\nType 1:<=\nType 2:>=\nType 3:=\n\n");
  int i,eqn_type[no_of_eqns],decision_variable_count=no_of_variables;
  double **B=(double**)malloc((no_of_eqns+1)*sizeof(double*));
  for(i=0;i<no_of_eqns+1;i++){
    B[i]=(double*)malloc((decision_variable_count)*sizeof(double));
  }
  double b[no_of_eqns+1];
  for(i=0;i<no_of_eqns;i++){
    printf("Enter the type of condition equation %d\n",i+1);
    scanf("%d",&eqn_type[i]);
    int j;
    if(eqn_type[i]==1||eqn_type[i]==3){
      printf("Enter the coefficients in condition equation %d\n",i+1);
      for(j=0;j<decision_variable_count;j++){
        scanf("%lf",&B[i][j]);
      }
      scanf("%lf",&b[i]);
    }else if(eqn_type[i]==2){
      no_of_variables+=1;
      printf("Enter the coefficients in condition equation %d\n",i+1);
      int j=0;
      for(;j<decision_variable_count;j++){
        scanf("%lf",&B[i][j]);
      }
      scanf("%lf",&b[i]);
    }
  }
  printf("Enter the coefficients of objective function\n");
  int j=0;
  for(j=0;j<decision_variable_count;j++){
    scanf("%lf",&B[no_of_eqns][j]);
  }
  scanf("%lf",&b[no_of_eqns]);
  double **T;
  int column[no_of_variables],row[no_of_eqns+10];
  T=(double**)malloc((no_of_eqns+11)*sizeof(double*));
  int counter=decision_variable_count;
  for(i=0;i<no_of_eqns;i++){
    T[i]=(double*)malloc((no_of_variables+1)*sizeof(double));
    int j;
    if(eqn_type[i]==1||eqn_type[i]==3){
      for(j=0;j<no_of_variables+1;j++){
        if(j<decision_variable_count){
          T[i][j]=B[i][j];
        }else if(j<no_of_variables){
          T[i][j]=0;
        }else{
          T[i][j]=b[i];
        }
      }
    }else if(eqn_type[i]==2){
      int flag=0;
      for(j=0;j<no_of_variables+1;j++){
        if(j<decision_variable_count){
          T[i][j]=B[i][j];
        }else if(j==counter){
          T[i][j]=-1;
          flag=1;
        }else if(j==no_of_variables){
          T[i][j]=b[i];
        }else{
          T[i][j]=0;
        }
      }
      if(flag==1)
        counter++;
    }
  }
  T[no_of_eqns]=(double*)malloc((no_of_variables+1)*sizeof(double));
  double sum[no_of_variables+1];
  for(j=0 ; j<no_of_variables+1 ; j++)
    sum[j] = 0;
  for(j=0;j<no_of_variables+1;j++){
    for(i=0;i<no_of_eqns;i++){
      if(eqn_type[i]!=1&&j<no_of_variables){
        sum[j]+=T[i][j];
      }else if(eqn_type[i]!=1&&j==no_of_variables){
        sum[j]+=T[i][j];
      }
    }
  }
  for(j=0;j<no_of_variables+1;j++){
    T[no_of_eqns][j]=-(sum[j]);
  }
  for(i=0;i<no_of_variables;i++){
    column[i]=10+i;
  }
  for(i=0;i<no_of_eqns;i++){
    row[i]=20+i;
  }
  int flag;
  printf("\nPhase-I Initial Table: \n");
  int k,l;
  for(l=0;l<no_of_variables;l++){
    printf("X%d\t\t",(column[l]%10)+1);
  }
  printf("1\n");
  for(k=0;k<no_of_eqns+1;k++){
    for(l=0;l<no_of_variables+1;l++){
        printf("%lf\t",T[k][l]);
    }
    if(k==no_of_eqns){
        printf("Z\n");
    }else{
        switch(row[k]/10){
      case 1:printf("X%d\n",(row[k]%10)+1);
        break;
      case 2:printf("Z%d\n",(row[k]%10)+1);
        break;
      }
    }
  }
  while(1){
    flag=0;
    int decider=simplex_method(T,no_of_eqns,no_of_variables,column,row,1,eqn_type);
    if(decider!=0){
      break;
    }
    for(i=0;i<no_of_variables;i++){
      if(fabs(T[no_of_eqns][i])<0.0001){
        T[no_of_eqns][i]=0;
      }
      if(T[no_of_eqns][i]<0){
        flag=1;
      }
    }
    if(flag==0){
      int flag1=0;
      for(i=0;i<no_of_eqns;i++){
        if(row[i]/10==2&&(eqn_type[row[i]%10]!=1)&&T[i][no_of_variables]!=0){
            flag1=1;
            printf("success\n");
        }
      }
      if(flag1==1){
        printf("\nThere is no solution for the following LPP.\n");
      }else if(flag1==0){
        print_result_phase_I(T,no_of_eqns,no_of_variables,column,row);
      }
      break;
    }
  }

    printf("\n\nPhase-II Initial table\n");
    double **temp=(double**)malloc((no_of_eqns+1)*sizeof(double*));
    int temp_column[no_of_variables];
    int count=0;
    for(i=0;i<=no_of_eqns;i++){
        temp[i]=(double*)malloc((no_of_variables+1)*sizeof(double));
    }
    for(j=0;j<no_of_variables;j++){
        if(column[j]/10==2&&eqn_type[column[j]%10]!=1){
            continue;
        }else{
            for(i=0;i<no_of_eqns;i++){
                temp[i][count]=T[i][j];
            }
            temp_column[count]=column[j];
            count++;
        }
    }
    for(i=0;i<no_of_eqns;i++){
        temp[i][count]=T[i][no_of_variables];
    }
    T=temp;
    for(i=0;i<no_of_variables;i++){
        column[i]=temp_column[i];
    }
  no_of_variables=count;
  int row_coefficient[no_of_eqns],column_coefficient[no_of_variables+1];
  for(i=0;i<no_of_eqns;i++){
        if(row[i]/10==1&&row[i]%10<decision_variable_count){
            row_coefficient[i]=B[no_of_eqns][row[i]%10];
        }else{
            row_coefficient[i]=0;
        }
    }

    for(i=0;i<no_of_variables;i++){
        if(column[i]/10==1&&column[i]%10<decision_variable_count){
            column_coefficient[i]=B[no_of_eqns][column[i]%10];
        }else{
            column_coefficient[i]=0;
        }
    }
    column_coefficient[no_of_variables]=b[no_of_eqns];

  for(j=0;j<=no_of_variables;j++){
        sum[j]=0;
    for(i=0;i<no_of_eqns;i++){
        sum[j]+=row_coefficient[i]*T[i][j];
    }
    if(j==no_of_variables){
        sum[j]+=column_coefficient[j];
    }else{
        sum[j]-=column_coefficient[j];
    }
    T[no_of_eqns][j]=sum[j];
  }

 for(l=0;l<no_of_variables;l++){
    printf("(%d)\t\t",column_coefficient[l]);
 }
 printf("\n");
  for(l=0;l<no_of_variables;l++){
    if(column[l]/10==1){
      printf("X%d\t\t",(column[l]%10)+1);
    }else if(column[l]/10==2){
      printf("Z%d\t\t",(column[l]%10)+1);
    }
  }
  printf("1\n");
  for(k=0;k<no_of_eqns+1;k++){
    for(l=0;l<no_of_variables+1;l++){
        printf("%lf\t",T[k][l]);
    }
    if(k==no_of_eqns){
        printf("Z\n");
    }else{
        switch(row[k]/10){
            case 1:printf("X%d",(row[k]%10)+1);
            break;
            case 2:printf("Z%d",(row[k]%10)+1);
            break;
        }
        printf("(%d\n)",row_coefficient[k]);
    }
  }
   while(1){
        flag=0;
        double min=100;
        int decider=simplex_method(T,no_of_eqns,no_of_variables,column,row,2,eqn_type);
        if(decider!=0){
            break;
        }
        for(i=0;i<no_of_variables;i++){
            if(fabs(T[no_of_eqns][i])<0.001){
                T[no_of_eqns][i]=0;
            }
            if(T[no_of_eqns][i]<0){
                flag=1;
            }
            if(T[no_of_eqns][i]<min){
                if(column[i]/10==2&&eqn_type[column[i]%10]!=1)
                    continue;
                min=T[no_of_eqns][i];
            }
        }
        if(min==0){
            print_result_phase_II(T,no_of_eqns,no_of_variables,column,row);
            flag=2;
            break;
        }
        if(flag==0){
            int flag1=0;
            for(i=0;i<no_of_eqns;i++){
                if(row[i]/10==2)
                    flag1=1;
            }
            if(flag1==1){
                printf("\nThere is no solution for the following LPP.\n");
            }else if(flag1==0){
                print_result_phase_II(T,no_of_eqns,no_of_variables,column,row);
            }
            break;
        }
    }

    printf("\n");
    if(flag==2){
        int decider=simplex_method(T,no_of_eqns,no_of_variables,column,row,2,eqn_type);
        if(decider==0){
            print_result_phase_II(T,no_of_eqns,no_of_variables,column,row);
        }
        printf("\nThere are infinite optimal solutions to this LPP");
    }

    while(1){


        double max=0,integer;
        int max_row=-1;
        for(i=0;i<no_of_eqns;i++){
            if((double)(T[i][no_of_variables]-floor(T[i][no_of_variables]))>(double)max){
                max=T[i][no_of_variables]<0?(T[i][no_of_variables]+1-(int)T[i][no_of_variables]):(T[i][no_of_variables]-(int)T[i][no_of_variables]);;
                max_row=i;
            }
        }
        double cc=modf((double)(max-floor(max)),&integer);
        bool tt=fabs(floor(cc)-cc)<0.00001L||fabs(ceil(cc)-cc)<0.00001L;
        if(tt==1){
          break;
        }else{
            row[no_of_eqns]=20+no_of_eqns;
            double temp1[no_of_variables+1];
            T[no_of_eqns+1]=(double*)malloc((no_of_variables+1)*sizeof(double));
            for(i=0;i<no_of_variables+1;i++){
                temp1[i]=T[no_of_eqns][i];
                T[no_of_eqns][i]=T[max_row][i]<0?-(T[max_row][i]+1-(int)T[max_row][i]):-(T[max_row][i]-(int)T[max_row][i]);
                T[no_of_eqns+1][i]=temp1[i];
            }
            no_of_eqns+=1;
        }


        printf("\nApplying Dual Simplex Method on updated LPP\n");
        printf("\nInitial Table: \n");
        for(l=0;l<no_of_variables;l++){
            printf("X%d\t\t",(column[l]%10)+1);
        }
        printf("Xb\n");
        for(k=0;k<no_of_eqns+1;k++){
          for(l=0;l<no_of_variables+1;l++){
            printf("%lf\t",T[k][l]);
          }
          if(k==no_of_eqns){
            printf("Z\n");
          }else{
              switch(row[k]/10){
                case 1:printf("X%d\n",(row[k]%10)+1);
                    break;
                case 2:printf("Z%d\n",(row[k]%10)+1);
                    break;
              }
          }
        }
        flag=0;

        for(i=0;i<no_of_variables;i++){
            if(T[no_of_eqns][i]<0)
                flag=1;
        }
        for(i=0;i<no_of_eqns;i++){
            if(T[i][no_of_variables]<0)
                flag=1;
        }
        if(flag==0){
            printf("The LPP does not satisfy conditions for Dual simplex method");
        }else if(flag==1){
            while(1){
                flag=0;
                double min=100;
                int decider=dual_simplex_method(T,no_of_eqns,no_of_variables,column,row);
                if(decider==1){
                    print_result(T,no_of_eqns,no_of_variables,column,row,no_of_variables);
                }
                if(decider!=0){
                    break;
                }
                for(i=0;i<no_of_variables;i++){
                    if(fabs(T[no_of_eqns][i])<0.001){
                        T[no_of_eqns][i]=0;
                    }
                    if(T[no_of_eqns][i]>0){
                        flag=1;
                    }
                    if(T[no_of_eqns][i]<min){
                        min=T[no_of_eqns][i];
                    }
                }
                if(min==0){
                    print_result(T,no_of_eqns,no_of_variables,column,row,no_of_variables);
                    flag=2;
                    break;
                }
                if(flag==0){
                    int flag1=0;
                    for(i=0;i<no_of_eqns;i++){
                        if(row[i]/10==2)
                            flag1=1;
                    }
                    if(flag1==1){
                        printf("\nThere is no solution for the following LPP.\n");
                    }else if(flag1==0){
                        print_result(T,no_of_eqns,no_of_variables,column,row,no_of_variables);
                    }
                    break;
                }
            }
        }
        printf("\n");
        if(flag==2){
            int decider=dual_simplex_method(T,no_of_eqns,no_of_variables,column,row);
            if(decider==0){
                print_result(T,no_of_eqns,no_of_variables,column,row,no_of_variables);
            }
            printf("\nThere are infinite optimal solutions to this LPP");
        }
    }

  return 1;
}

