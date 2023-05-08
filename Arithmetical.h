/* Kelei999999(WangLiang) all rights reserved.  You may use this software
 * and any derivatives exclusively with Kelei999999(WangLiang) products.
 *
 * THIS SOFTWARE IS SUPPLIED BY Kelei999999(WangLiang) "AS IS".  NO WARRANTIES, WHETHER
 * EXPRESS, IMPLIED OR STATUTORY, APPLY TO THIS SOFTWARE, INCLUDING ANY IMPLIED
 * WARRANTIES OF NON-INFRINGEMENT, MERCHANTABILITY, AND FITNESS FOR A
 * PARTICULAR PURPOSE, OR ITS INTERACTION WITH Kelei999999(WangLiang) PRODUCTS, COMBINATION
 * WITH ANY OTHER PRODUCTS, OR USE IN ANY APPLICATION.
 *
 * IN NO EVENT WILL Kelei999999(WangLiang) BE LIABLE FOR ANY INDIRECT, SPECIAL, PUNITIVE,
 * INCIDENTAL OR CONSEQUENTIAL LOSS, DAMAGE, COST OR EXPENSE OF ANY KIND
 * WHATSOEVER RELATED TO THE SOFTWARE, HOWEVER CAUSED, EVEN IF Kelei999999(WangLiang) HAS
 * BEEN ADVISED OF THE POSSIBILITY OR THE DAMAGES ARE FORESEEABLE.  TO THE
 * FULLEST EXTENT ALLOWED BY LAW, Kelei999999(WangLiang)'S TOTAL LIABILITY ON ALL CLAIMS
 * IN ANY WAY RELATED TO THIS SOFTWARE WILL NOT EXCEED THE AMOUNT OF FEES, IF
 * ANY, THAT YOU HAVE PAID DIRECTLY TO Kelei999999(WangLiang) FOR THIS SOFTWARE.
 *
 * Kelei999999(WangLiang) PROVIDES THIS SOFTWARE CONDITIONALLY UPON YOUR ACCEPTANCE OF THESE
 * TERMS.
 */

/*
 * File: Arithmetical.h
 * Author: Kelei999999(WangLiang)
 * Created on: 2014年4月12日
 * Revision history: 1.1
 */

/*
 * This is a guard condition so that contents of this file are not included
 * more than once.
 */

#ifndef ARITHMETICAL_H_
#define ARITHMETICAL_H_
#include "Typedef.h"
typedef struct
{
	volatile real32 *KP;					/*本次PID运算的P值*/
	volatile real32 *KI;					/*本次PID运算的I值*/
	volatile real32 *current_Value;			/*本次PID运算给定量的当前值*/
	volatile real32 *target_Value;			/*本次PID运算给定量的目标值*/
	volatile real32 *I_Value;				/*上次PID运算的积分值（调用者需要自行保存好）*/
	volatile real32 *ReferenceRatio;		/*本次PID运算时误差的比例系数，其值最好是给定量最大值的50%的倒数*/
}pid;
/*
 * 增量式PID运算结构体
 */
typedef struct
{
	volatile real32 *KP;					/*本次PID运算的P值*/
	volatile real32 *KI;					/*本次PID运算的I值*/
	volatile real32 *KD;					/*本次PID运算的D值*/
	volatile real32 *current_Value;			/*本次PID运算给定量的当前值*/
	volatile real32 *target_Value;			/*本次PID运算给定量的目标值*/
	volatile real32 *error_Last;			/*PID运算的前次误差（调用者需要自行保存好）*/
	volatile real32 *Delta_Last;			/*PID运算的误差增量的导数（调用者需要自行保存好）*/
	volatile real32 *PI_OUT;				/*PID的输出值（已经被累加了）*/
	volatile real32 *Regulator_out;			/*调节器最终的输出值（在用户使用中可能已经被限幅了，在没有限幅之前其值=PI_OUT），用于防止积分饱和*/
	volatile real32 *ReferenceRatio;		/*PID运算时误差的比例系数，其值最好是给定量最大值的50%的倒数*/
}pidi;
/*
 * 增量式PID运算结构体
 */
typedef struct
{
	volatile real32 *KP;					/*本次PID运算的P值*/
	volatile real32 *KI;					/*本次PID运算的I值*/
	volatile real32 *KD;					/*本次PID运算的D值*/
	volatile real32 *current_Value;			/*本次PID运算给定量的当前值*/
	volatile real32 *target_Value;			/*本次PID运算给定量的目标值*/
	volatile real32 *error_Last;			/*PID运算的前次误差（调用者需要自行保存好）*/
	volatile real32 *Delta_Last;			/*PID运算的误差增量的导数（调用者需要自行保存好）*/
	volatile real32 *PI_OUT;				/*PID的输出值（已经被累加了）*/
	volatile real32 *Regulator_out;			/*调节器最终的输出值（在用户使用中可能已经被限幅了，在没有限幅之前其值=PI_OUT），用于防止积分饱和*/
	volatile real32 *ReferenceRatio;		/*PID运算时误差的比例系数，其值最好是给定量最大值的50%的倒数*/
	volatile real32 *Delta_Cur;				/*PID运算时的电流差分值*/
	volatile real32 *KC;					/*本次PID运算的C值*/
}pidic;
/*
 * 模糊控制增量式PID运算结构体
 */
typedef struct
{
	volatile real32 *KP;					/*本次PID运算的P值*/
	volatile real32 *KI;					/*本次PID运算的I值*/
	volatile real32 *KD;					/*本次PID运算的D值*/
	volatile real32 *current_Value;			/*本次PID运算给定量的当前值*/
	volatile real32 *target_Value;			/*本次PID运算给定量的目标值*/
	volatile real32 *error_Last;			/*PID运算的前次误差（调用者需要自行保存好）*/
	volatile real32 *Delta_Last;			/*PID运算的误差增量的导数（调用者需要自行保存好）*/
	volatile real32 *PI_OUT;				/*PID的输出值（已经被累加了）*/
	volatile real32 *Regulator_out;			/*调节器最终的输出值（在用户使用中可能已经被限幅了，在没有限幅之前其值=PI_OUT），用于防止积分饱和*/
	volatile real32 *ReferenceRatio;		/*PID运算时误差的比例系数，其值最好是给定量最大值的50%的倒数*/
	volatile real32 *Delta_Cur;				/*PID运算时的电流差分值*/
	volatile real32 *KC;					/*本次PID运算的C值*/
	volatile real32 *KFuzz;					/*本次PID运算的模糊禁带宽度*/
}pidif;
/*
 * 模糊控制增量式PID V1运算结构体
 */
typedef struct
{
	volatile real32 *KP;					/*本次PID运算的P值*/
	volatile real32 *KI;					/*本次PID运算的I值*/
	volatile real32 *KD;					/*本次PID运算的D值*/
	volatile real32 *current_Value;			/*本次PID运算给定量的当前值*/
	volatile real32 *target_Value;			/*本次PID运算给定量的目标值*/
	volatile real32 *error_Last;			/*PID运算的前次误差（调用者需要自行保存好）*/
	volatile real32 *Delta_Last;			/*PID运算的误差增量的导数（调用者需要自行保存好）*/
	volatile real32 *PI_OUT;				/*PID的输出值（已经被累加了）*/
	volatile real32 *ReferenceRatio;		/*PID运算时误差的比例系数，其值最好是给定量最大值的50%的倒数*/
	volatile real32 *Delta_Cur;				/*PID运算时的电流差分值*/
	volatile real32 *KC;					/*本次PID运算的C值*/
	volatile real32 *KFuzz;					/*本次PID运算的模糊禁带宽度*/
	volatile real32 UpperLimit;				/*PID输出上限值*/
	volatile real32 LowerLimit;				/*PID输出下限值*/
}pidifv1;
/*
 * 模糊控制增量式PID V2运算结构体
 */
typedef struct
{
	volatile real32 *KP;					/*本次PID运算的P值*/
	volatile real32 *KI;					/*本次PID运算的I值*/
	volatile real32 *KD;					/*本次PID运算的D值*/
	volatile real32 *current_Value;			/*本次PID运算给定量的当前值*/
	volatile real32 *target_Value;			/*本次PID运算给定量的目标值*/
	volatile real32 *error_Last;			/*PID运算的前次误差（调用者需要自行保存好）*/
	volatile real32 *Delta_Last;			/*PID运算的误差增量的导数（调用者需要自行保存好）*/
	volatile real32 *PI_OUT;				/*PID的输出值（已经被累加了）*/
	volatile real32 *ReferenceRatio;		/*PID运算时误差的比例系数，其值最好是给定量最大值的50%的倒数*/
	volatile real32 *Delta_Cur;				/*PID运算时的电流差分值*/
	volatile real32 *KC;					/*本次PID运算的C值*/
	volatile real32 *KFuzz;					/*本次PID运算的模糊禁带宽度*/
	volatile real32 UpperLimit;				/*PID输出上限值*/
	volatile real32 LowerLimit;				/*PID输出下限值*/
	volatile real32 ErrorLimit;				/*PID误差限值正负范围*/
}pidifv2;
/*
 * 模糊控制增量式2p2z PID运算结构体
 */
typedef struct
{
	struct
	{
		volatile real32 *target_Value;			/*本次PID运算给定量的目标值*/
		volatile real32 *current_Value;			/*本次PID运算给定量的当前值*/
		volatile real32 *Delta_Cur;				/*PID运算时的电流差分值，配合KC使用，不使用时应置零*/
		volatile real32 *KP_1p;					/*本次PID运算的P值*/
		volatile real32 *KI_1p;					/*本次PID运算的I值*/
		volatile real32 *KD_1p;					/*本次PID运算的D值*/
		volatile real32 *KC_1p;					/*本次PID运算的C值(电流补偿值，仅在电流PI中有效)*/
		volatile real32 *KF_1p;					/*本次PID运算的模糊禁带宽度*/
		volatile real32 *KP_2p;					/*本次PID运算的P值*/
		volatile real32 *KI_2p;					/*本次PID运算的I值*/
		volatile real32 *KD_2p;					/*本次PID运算的D值*/
		volatile real32 *KC_2p;					/*本次PID运算的C值(电流补偿值，仅在电流PI中有效)*/
		volatile real32 *KF_2p;					/*本次PID运算的模糊禁带宽度*/
		volatile real32 RefRatio;				/*PID运算时误差的比例系数，其值最好是给定量最大值的50%的倒数*/
	}input;
	struct
	{
		volatile real32 error_Last_1p;			/*PID运算的前次误差（调用者需要自行保存好）*/
		volatile real32 Delta_Last_1p;			/*PID运算的误差增量的导数（调用者需要自行保存好）*/
		volatile real32 error_Last_2p;			/*PID运算的前次误差（调用者需要自行保存好）*/
		volatile real32 Delta_Last_2p;			/*PID运算的误差增量的导数（调用者需要自行保存好）*/
		volatile real32 PI_OUT;					/*PID的输出值（已经被累加了）*/
		volatile real32 Regulator_out;			/*调节器最终的输出值（在用户使用中可能已经被限幅了，在没有限幅之前其值=PI_OUT），用于防止积分饱和*/
	}record;
}pidif_2p2z;
/*
 * 均方根值计算结构体
 */
typedef struct
{
	struct 									/*输入值，以下数据由用户提供*/
	{
		volatile real32 *INS;				/*需要RMS计算信号的瞬时值*/
		volatile real32 Threshold;			/*需要RMS计算信号的零位门限阈值[-Threshold,Threshold]*/
		volatile uint32 Count_MAX;			/*RMS计算最大允许累加次数*/
	}input;
	struct 									/*用于RMS计算的中间量,以下数据用户不允许改变*/
	{
		volatile real32 RMS_Value;			/*RMS计算的结果*/
		volatile real32 RMS_ACC;			/*RMS计算累加值*/
		volatile uint8 ZERO;				/*需要RMS计算信号的零位标志*/
		volatile uint8 LOCK;				/*需要RMS计算信号的零位锁定标志*/
		volatile uint32 TIME;				/*RMS计算累加次数*/
	}temp;
}rms_process;
//****************************************************************************
// @函数名                   complex CLARKE(real32_T A,real32_T B,real32_T C)
//----------------------------------------------------------------------------
// @描述                       函数实现了ckarke等幅变换
//
//----------------------------------------------------------------------------
// @输入                      三相静止坐标下的瞬时值
//
//----------------------------------------------------------------------------
// @输出                       两相静止坐标下的瞬时值
//
//----------------------------------------------------------------------------
// @返回值		complex
//			complex定义为:
//			typedef struct
//			{
//				float alafa;
//				float beita;
//			}complex;
//
//----------------------------------------------------------------------------
// @日期                       2017年3月21日
//
//****************************************************************************
extern complex CLARKE(real32_T A,real32_T B,real32_T C);
//****************************************************************************
// @函数名			real32 PID(pid *Pid)
//----------------------------------------------------------------------------
// @描述			函数实现了PID运算
//					误差ERROR = 目标值 - 当前值
//----------------------------------------------------------------------------
// @输入			pid参数，详见pid结构体
//
//----------------------------------------------------------------------------
// @输出			PID运算后的值
//
//----------------------------------------------------------------------------
// @返回值			real32 PID运算后的值
//----------------------------------------------------------------------------
// @日期			2017年4月20日
//
//****************************************************************************
extern real32 PID(pid *Pid);
//****************************************************************************
// @函数名			real32 PIDI(pidi *Pid)
//----------------------------------------------------------------------------
// @描述			函数实现了增量式PID运算
//					误差ERROR = 目标值 - 当前值
//----------------------------------------------------------------------------
// @输入			pid参数，详见pidi结构体
//
//----------------------------------------------------------------------------
// @输出			PID运算后的值
//
//----------------------------------------------------------------------------
// @返回值			real32 PID运算后的值
//----------------------------------------------------------------------------
// @日期			2017年4月20日
//
//****************************************************************************
extern real32 PIDI(pidi *Pid);
//****************************************************************************
// @函数名			real32 PIDC(pidc *Pid)
//----------------------------------------------------------------------------
// @描述			函数实现了带电流差分补偿的增量式PID运算
//					误差ERROR = 目标值 - 当前值
//----------------------------------------------------------------------------
// @输入			pid参数，详见pidc结构体
//
//----------------------------------------------------------------------------
// @输出			PID运算后的值
//
//----------------------------------------------------------------------------
// @返回值			real32 PID运算后的值
//----------------------------------------------------------------------------
// @日期			2017年4月20日
//
//****************************************************************************
extern real32 PIDIC(pidic *Pid);
//****************************************************************************
// @函数名			real32 PIDICF(pidif *Pid)
//----------------------------------------------------------------------------
// @描述			函数实现了带电流差分补偿和模糊禁带宽度的增量式PID运算
//					当给定值小于2倍禁带宽度时，禁带宽度将不起作用
//					误差ERROR = 目标值 - 当前值
//----------------------------------------------------------------------------
// @输入			pid参数，详见pidif结构体
//
//----------------------------------------------------------------------------
// @输出			PID运算后的值
//
//----------------------------------------------------------------------------
// @返回值			real32 PID运算后的值
//----------------------------------------------------------------------------
// @日期			2017年4月20日
//
//****************************************************************************
extern real32 PIDICF(pidif *Pid);
//****************************************************************************
// @函数名			real32 PIDICF_INV(pidif *Pid)
//----------------------------------------------------------------------------
// @描述			函数实现了带电流差分补偿和模糊禁带宽度的增量式反向PID运算
//					当给定值小于2倍禁带宽度时，禁带宽度将不起作用
//					误差ERROR = 当前值 - 目标值
//----------------------------------------------------------------------------
// @输入			pid参数，详见pidif结构体
//
//----------------------------------------------------------------------------
// @输出			PID运算后的值
//
//----------------------------------------------------------------------------
// @返回值			real32 PID运算后的值
//----------------------------------------------------------------------------
// @日期			2017年4月20日
//
//****************************************************************************
extern real32 PIDICF_INV(pidif *Pid);
//****************************************************************************
// @函数名			real32 SYNC_RMS_V1(rms_process *RMS)
//----------------------------------------------------------------------------
// @描述			函数实现了伪同步RMS计算
//
//----------------------------------------------------------------------------
// @输入			rms_process *RMS:详见rms_process结构体
//
//----------------------------------------------------------------------------
// @输出			无
//
//----------------------------------------------------------------------------
// @返回值			RMS的值
//----------------------------------------------------------------------------
// @日期			2017年9月13日
//
//****************************************************************************
extern real32 SYNC_RMS_V1(rms_process *RMS);
//****************************************************************************
// @函数名			real32 PIDICFV1(pidifv1 *Pid)
//----------------------------------------------------------------------------
// @描述			函数实现了带电流差分补偿和模糊禁带宽度的增量式PID运算
//					当给定值小于2倍禁带宽度时，禁带宽度将不起作用
//					误差ERROR = 目标值 - 当前值
//----------------------------------------------------------------------------
// @输入			pid参数，详见pidifv1结构体
//
//----------------------------------------------------------------------------
// @输出			PID运算后的值
//
//----------------------------------------------------------------------------
// @返回值			real32 PID运算后的值
//----------------------------------------------------------------------------
// @日期			2017年4月20日
//
//****************************************************************************
extern real32 PIDICFV1(pidifv1 *Pid);
//****************************************************************************
// @函数名			real32 PIDICF_INV1(pidifv1 *Pid)
//----------------------------------------------------------------------------
// @描述			函数实现了带电流差分补偿和模糊禁带宽度的增量式反向PID运算
//					当给定值小于2倍禁带宽度时，禁带宽度将不起作用
//					误差ERROR = 当前值 - 目标值
//----------------------------------------------------------------------------
// @输入			pid参数，详见pidifv1结构体
//
//----------------------------------------------------------------------------
// @输出			PID运算后的值
//
//----------------------------------------------------------------------------
// @返回值			real32 PID运算后的值
//----------------------------------------------------------------------------
// @日期			2017年4月26日
//
//****************************************************************************
extern real32 PIDICF_INV1(pidifv1 *Pid);
//****************************************************************************
// @函数名			real32 PIDICFV2(pidifv2 *Pid)
//----------------------------------------------------------------------------
// @描述			函数实现了带电流差分补偿和模糊禁带宽度的增量式PID运算
//					当给定值小于2倍禁带宽度时，禁带宽度将不起作用
//					误差ERROR = 目标值 - 当前值
//----------------------------------------------------------------------------
// @输入			pid参数，详见pidifv2结构体
//
//----------------------------------------------------------------------------
// @输出			PID运算后的值
//
//----------------------------------------------------------------------------
// @返回值			real32 PID运算后的值
//----------------------------------------------------------------------------
// @日期			2017年4月20日
//
//****************************************************************************
extern real32 PIDICFV2(pidifv2 *Pid);
//****************************************************************************
// @函数名			real32 PIDICF_INV2(pidifv2 *Pid)
//----------------------------------------------------------------------------
// @描述			函数实现了带电流差分补偿和模糊禁带宽度的增量式反向PID运算
//					当给定值小于2倍禁带宽度时，禁带宽度将不起作用
//					误差ERROR = 当前值 - 目标值
//----------------------------------------------------------------------------
// @输入			pid参数，详见pidifv2结构体
//
//----------------------------------------------------------------------------
// @输出			PID运算后的值
//
//----------------------------------------------------------------------------
// @返回值			real32 PID运算后的值
//----------------------------------------------------------------------------
// @日期			2017年4月26日
//
//****************************************************************************
extern real32 PIDICF_INV2(pidifv2 *Pid);
#endif /* ARITHMETICAL_H_ */
