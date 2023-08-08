/* Hunter.ORG all rights reserved.  You may use this software
 * and any derivatives exclusively with Hunter.ORG products.
 *
 * THIS SOFTWARE IS SUPPLIED BY Hunter.ORG "AS IS".  NO WARRANTIES, WHETHER
 * EXPRESS, IMPLIED OR STATUTORY, APPLY TO THIS SOFTWARE, INCLUDING ANY IMPLIED
 * WARRANTIES OF NON-INFRINGEMENT, MERCHANTABILITY, AND FITNESS FOR A
 * PARTICULAR PURPOSE, OR ITS INTERACTION WITH Hunter.ORG PRODUCTS, COMBINATION
 * WITH ANY OTHER PRODUCTS, OR USE IN ANY APPLICATION.
 *
 * IN NO EVENT WILL Hunter.ORG BE LIABLE FOR ANY INDIRECT, SPECIAL, PUNITIVE,
 * INCIDENTAL OR CONSEQUENTIAL LOSS, DAMAGE, COST OR EXPENSE OF ANY KIND
 * WHATSOEVER RELATED TO THE SOFTWARE, HOWEVER CAUSED, EVEN IF Hunter.ORG HAS
 * BEEN ADVISED OF THE POSSIBILITY OR THE DAMAGES ARE FORESEEABLE.  TO THE
 * FULLEST EXTENT ALLOWED BY LAW, Hunter.ORG'S TOTAL LIABILITY ON ALL CLAIMS
 * IN ANY WAY RELATED TO THIS SOFTWARE WILL NOT EXCEED THE AMOUNT OF FEES, IF
 * ANY, THAT YOU HAVE PAID DIRECTLY TO Hunter.ORG FOR THIS SOFTWARE.
 *
 * Hunter.ORG PROVIDES THIS SOFTWARE CONDITIONALLY UPON YOUR ACCEPTANCE OF THESE
 * TERMS.
 */

/*
 * File: Hunter_OS_Arith.c
 * Author: Hunter.ORG
 * Created on: 2014年4月12日
 * Revision history: 1.1
 */

/*
 * This is a guard condition so that contents of this file are not included
 * more than once.
 */
#include "Hunter_OS_Arith.h"
#include <math.h>
//****************************************************************************
// @函数名			complex CLARKE(real32_T A,real32_T B,real32_T C)
//----------------------------------------------------------------------------
// @描述			函数实现了ckarke等幅变换
//
//----------------------------------------------------------------------------
// @输入			三相静止坐标下的瞬时值
//
//----------------------------------------------------------------------------
// @输出			两相静止坐标下的瞬时值
//
//----------------------------------------------------------------------------
// @返回值			complex
//					complex定义为:
//						typedef struct
//						{
//							float alafa;
//							float beita;
//						}complex;
//
//----------------------------------------------------------------------------
// @日期			2017年3月21日
//
//****************************************************************************
static complex CLARKE(real32 A,real32 B,real32 C)
{
	complex Temp;
	Temp.alpha = (A - (B + C) * 0.5f) * 0.666666666666f;//0.666666666666f = (2/3)
	Temp.beta = ((B - C) * 0.8660254f) * 0.666666666666f;
	return Temp;
}
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
static real32 PID(pid *Pid)
{
	static real32	VF32_TEM_PI_ERROR;
	VF32_TEM_PI_ERROR = (*Pid->target_Value - *Pid->current_Value) * *Pid->ReferenceRatio;
	if(VF32_TEM_PI_ERROR > 1.1f)
		VF32_TEM_PI_ERROR = 1.1f;
	else if(VF32_TEM_PI_ERROR < -1.1f)
		VF32_TEM_PI_ERROR = -1.1f;

	*Pid->I_Value += VF32_TEM_PI_ERROR * 0.0005f;

	if(VF32_TEM_PI_ERROR > 0.2f)
		VF32_TEM_PI_ERROR = 0.2f;
	else if(VF32_TEM_PI_ERROR < -0.2f)
		VF32_TEM_PI_ERROR = -0.2f;

	if(*Pid->I_Value > 1.0f)
		*Pid->I_Value = 1.0f;
	else if(*Pid->I_Value < -1.0f)
		*Pid->I_Value = -1.0f;

	return (*Pid->KP * VF32_TEM_PI_ERROR + *Pid->KI * *Pid->I_Value);
}

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
//extern real32 					VF32_TEM_TEMP1,VF32_TEM_TEMP2,VF32_TEM_TEMP3;
static real32 PIDI(pidi *Pid)
{
	real32	VF32_TEM_PI_ERROR;				/*PID运算的当前误差值*/
	real32	VF32_TEM_PI_ERROR_Delta;		/*PID运算的误差值增量*/
	real32	VF32_TEM_PI_pERROR_Delta;		/*PID运算的误差值增量的导数*/
	real32	PIDI_TEMP;
	//--------------------------------------获取当前误差
	VF32_TEM_PI_ERROR = (*Pid->target_Value - *Pid->current_Value) * *Pid->ReferenceRatio;
	//--------------------------------------PID误差限幅
	if(VF32_TEM_PI_ERROR > 0.2f)
		VF32_TEM_PI_ERROR = 0.2f;
	else if(VF32_TEM_PI_ERROR < -0.2f)
		VF32_TEM_PI_ERROR = -0.2f;
	//--------------------------------------获取当前误差的增量
	VF32_TEM_PI_ERROR_Delta = VF32_TEM_PI_ERROR - *Pid->error_Last;
	//--------------------------------------获取当前误差增量的导数
	VF32_TEM_PI_pERROR_Delta = VF32_TEM_PI_ERROR_Delta - *Pid->Delta_Last;
	//--------------------------------------采用梯形积分法求取积分值
	PIDI_TEMP = *Pid->KI * (VF32_TEM_PI_ERROR + *Pid->error_Last) - (*Pid->PI_OUT - *Pid->Regulator_out);
	//--------------------------------------PID运算
	*Pid->PI_OUT += (*Pid->KP * VF32_TEM_PI_ERROR_Delta + PIDI_TEMP + *Pid->KD * VF32_TEM_PI_pERROR_Delta);
	//--------------------------------------存储当前的误差和上次的误差
	*Pid->Delta_Last = VF32_TEM_PI_ERROR_Delta;
	*Pid->error_Last = VF32_TEM_PI_ERROR;
	//--------------------------------------PID输出
	return *Pid->PI_OUT;
}
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
//extern real32 					VF32_TEM_TEMP1,VF32_TEM_TEMP2,VF32_TEM_TEMP3;
static real32 PIDIC(pidic *Pid)
{
	real32	VF32_TEM_PI_ERROR;				/*PID运算的当前误差值*/
	real32	VF32_TEM_PI_ERROR_Delta;		/*PID运算的误差值增量*/
	real32	VF32_TEM_PI_pERROR_Delta;		/*PID运算的误差值增量的导数*/
	real32	PIDI_TEMP;
	//--------------------------------------获取当前误差
	VF32_TEM_PI_ERROR = (*Pid->target_Value - *Pid->current_Value) * *Pid->ReferenceRatio;
	//--------------------------------------PID误差限幅
	if(VF32_TEM_PI_ERROR > 0.8f)
		VF32_TEM_PI_ERROR = 0.8f;
	else if(VF32_TEM_PI_ERROR < -0.8f)
		VF32_TEM_PI_ERROR = -0.8f;
	//--------------------------------------获取当前误差的增量
	VF32_TEM_PI_ERROR_Delta = VF32_TEM_PI_ERROR - *Pid->error_Last;
	//--------------------------------------获取当前误差增量的导数
	VF32_TEM_PI_pERROR_Delta = VF32_TEM_PI_ERROR_Delta - *Pid->Delta_Last;
	//--------------------------------------采用梯形积分法求取积分值
	PIDI_TEMP = *Pid->KI * (VF32_TEM_PI_ERROR + *Pid->error_Last) - (*Pid->PI_OUT - *Pid->Regulator_out);
	//--------------------------------------PID运算
	*Pid->PI_OUT += (*Pid->KP * VF32_TEM_PI_ERROR_Delta + PIDI_TEMP + *Pid->KD * VF32_TEM_PI_pERROR_Delta + *Pid->Delta_Cur * *Pid->KC);
	//--------------------------------------存储当前的误差和上次的误差
	*Pid->Delta_Last = VF32_TEM_PI_ERROR_Delta;
	*Pid->error_Last = VF32_TEM_PI_ERROR;
	//--------------------------------------PID输出
	return *Pid->PI_OUT;
}
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
static real32 PIDICF(pidif *Pid)
{
	real32	VF32_TEM_PI_ERROR;				/*PID运算的当前误差值*/
	real32	VF32_TEM_PI_ERROR_Delta;		/*PID运算的误差值增量*/
	real32	VF32_TEM_PI_pERROR_Delta;		/*PID运算的误差值增量的导数*/
	real32	VF32_TEM_TEMP;
	//--------------------------------------获取当前误差
	VF32_TEM_TEMP = (*Pid->target_Value - *Pid->current_Value);
	//--------------------------------------禁带计算
	if(*Pid->target_Value < *Pid->KFuzz * 1.2f)//当给定值小于1.2f倍禁带宽度时，不对禁带宽度进行干预
	{
		VF32_TEM_PI_ERROR = VF32_TEM_TEMP * *Pid->ReferenceRatio;
	}
	else
	{
		if(VF32_TEM_TEMP >= *Pid->KFuzz)
			VF32_TEM_PI_ERROR = (VF32_TEM_TEMP - *Pid->KFuzz) * *Pid->ReferenceRatio;
		else if(VF32_TEM_TEMP <= -(*Pid->KFuzz))
			VF32_TEM_PI_ERROR = (VF32_TEM_TEMP + *Pid->KFuzz) * *Pid->ReferenceRatio;
		else
			VF32_TEM_PI_ERROR = 0;
	}
	//--------------------------------------PID误差限幅
	if(VF32_TEM_PI_ERROR > 0.4f)
		VF32_TEM_PI_ERROR = 0.4f;
	else if(VF32_TEM_PI_ERROR < -0.4f)
		VF32_TEM_PI_ERROR = -0.4f;
	//--------------------------------------获取当前误差的增量
	VF32_TEM_PI_ERROR_Delta = VF32_TEM_PI_ERROR - *Pid->error_Last;
	//--------------------------------------获取当前误差增量的导数
	VF32_TEM_PI_pERROR_Delta = VF32_TEM_PI_ERROR_Delta - *Pid->Delta_Last;
	//--------------------------------------采用梯形积分法求取积分值
	VF32_TEM_TEMP = *Pid->KI * (VF32_TEM_PI_ERROR + *Pid->error_Last) - (*Pid->PI_OUT - *Pid->Regulator_out);
	//--------------------------------------PID运算
	*Pid->PI_OUT += (*Pid->KP * VF32_TEM_PI_ERROR_Delta + VF32_TEM_TEMP + *Pid->KD * VF32_TEM_PI_pERROR_Delta + *Pid->Delta_Cur * *Pid->KC);
	//--------------------------------------存储当前的误差和上次的误差
	*Pid->Delta_Last = VF32_TEM_PI_ERROR_Delta;
	*Pid->error_Last = VF32_TEM_PI_ERROR;
	//--------------------------------------PID输出
	return *Pid->PI_OUT;
}
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
// @日期			2017年4月26日
//
//****************************************************************************
static real32 PIDICF_INV(pidif *Pid)
{
	real32	VF32_TEM_PI_ERROR;				/*PID运算的当前误差值*/
	real32	VF32_TEM_PI_ERROR_Delta;		/*PID运算的误差值增量*/
	real32	VF32_TEM_PI_pERROR_Delta;		/*PID运算的误差值增量的导数*/
	real32	VF32_TEM_TEMP;
	//--------------------------------------获取当前误差
	VF32_TEM_TEMP = (*Pid->current_Value - *Pid->target_Value);
	//--------------------------------------禁带计算
	if(*Pid->target_Value < *Pid->KFuzz * 1.2f)//当给定值小于1.2f倍禁带宽度时，不对禁带宽度进行干预
	{
		VF32_TEM_PI_ERROR = VF32_TEM_TEMP * *Pid->ReferenceRatio;
	}
	else
	{
		if(VF32_TEM_TEMP >= *Pid->KFuzz)
			VF32_TEM_PI_ERROR = (VF32_TEM_TEMP - *Pid->KFuzz) * *Pid->ReferenceRatio;
		else if(VF32_TEM_TEMP <= -(*Pid->KFuzz))
			VF32_TEM_PI_ERROR = (VF32_TEM_TEMP + *Pid->KFuzz) * *Pid->ReferenceRatio;
		else
			VF32_TEM_PI_ERROR = 0;
	}
	//--------------------------------------PID误差限幅
	if(VF32_TEM_PI_ERROR > 0.4f)
		VF32_TEM_PI_ERROR = 0.4f;
	else if(VF32_TEM_PI_ERROR < -0.4f)
		VF32_TEM_PI_ERROR = -0.4f;
	//--------------------------------------获取当前误差的增量
	VF32_TEM_PI_ERROR_Delta = VF32_TEM_PI_ERROR - *Pid->error_Last;
	//--------------------------------------获取当前误差增量的导数
	VF32_TEM_PI_pERROR_Delta = VF32_TEM_PI_ERROR_Delta - *Pid->Delta_Last;
	//--------------------------------------采用梯形积分法求取积分值
	VF32_TEM_TEMP = *Pid->KI * (VF32_TEM_PI_ERROR + *Pid->error_Last) - (*Pid->PI_OUT - *Pid->Regulator_out);
	//--------------------------------------PID运算
	*Pid->PI_OUT += (*Pid->KP * VF32_TEM_PI_ERROR_Delta + VF32_TEM_TEMP + *Pid->KD * VF32_TEM_PI_pERROR_Delta + *Pid->Delta_Cur * *Pid->KC);
	//--------------------------------------存储当前的误差和上次的误差
	*Pid->Delta_Last = VF32_TEM_PI_ERROR_Delta;
	*Pid->error_Last = VF32_TEM_PI_ERROR;
	//--------------------------------------PID输出
	return *Pid->PI_OUT;
}
//****************************************************************************
// @函数名		real32 SYNC_RMS_V1(rms_process *RMS)
//----------------------------------------------------------------------------
// @描述                       函数实现了伪同步RMS计算
//
//----------------------------------------------------------------------------
// @输入			rms_process *RMS:详见rms_process结构体
//
//----------------------------------------------------------------------------
// @输出
//
//----------------------------------------------------------------------------
// @返回值		RMS的值
//----------------------------------------------------------------------------
// @日期                       2017年5月28日
//
//****************************************************************************
static real32 SYNC_RMS_V1(rms_process *RMS)
{
	real32 temp = *RMS->input.INS;
	//--------------------------------------------找过零点
	if(temp > RMS->input.Threshold)
	{
		if(0 == RMS->temp.LOCK)
		{
			RMS->temp.ZERO = 1;
			RMS->temp.LOCK = 1;
		}
	}
	else if(temp < -RMS->input.Threshold)
	{
		if(1 == RMS->temp.LOCK)
		{
			RMS->temp.LOCK = 0;
		}
	}
	//-------------------------------------------计算平方累和
	RMS->temp.RMS_ACC += (temp * temp);
	//-------------------------------------------记录累和次数
	RMS->temp.TIME += 1;
	//-------------------------------------------
	if(RMS->temp.TIME >= RMS->input.Count_MAX)
		RMS->temp.ZERO |= 1;
	//-------------------------------------------计算均方根
	if(1 == RMS->temp.ZERO)
	{
		RMS->temp.ZERO = 0;
		if(RMS->temp.TIME > 0)
			RMS->temp.RMS_Value = sqrtf(RMS->temp.RMS_ACC / RMS->temp.TIME);
		else
			RMS->temp.RMS_Value = 0;
		RMS->temp.TIME = 0;
		RMS->temp.RMS_ACC = 0;
	}
	return RMS->temp.RMS_Value;
}
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
static real32 PIDICFV1(pidifv1 *Pid)
{
	real32	VF32_TEM_PI_ERROR;				/*PID运算的当前误差值*/
	real32	VF32_TEM_PI_ERROR_Delta;		/*PID运算的误差值增量*/
	real32	VF32_TEM_PI_pERROR_Delta;		/*PID运算的误差值增量的导数*/
	real32	VF32_TEM_TEMP;
	//--------------------------------------获取当前误差
	VF32_TEM_TEMP = (*Pid->target_Value - *Pid->current_Value);
	//--------------------------------------禁带计算
	if(*Pid->target_Value < *Pid->KFuzz * 1.2f)//当给定值小于1.2f倍禁带宽度时，不对禁带宽度进行干预
	{
		VF32_TEM_PI_ERROR = VF32_TEM_TEMP * *Pid->ReferenceRatio;
	}
	else
	{
		if(VF32_TEM_TEMP >= *Pid->KFuzz)
			VF32_TEM_PI_ERROR = (VF32_TEM_TEMP - *Pid->KFuzz) * *Pid->ReferenceRatio;
		else if(VF32_TEM_TEMP <= -(*Pid->KFuzz))
			VF32_TEM_PI_ERROR = (VF32_TEM_TEMP + *Pid->KFuzz) * *Pid->ReferenceRatio;
		else
			VF32_TEM_PI_ERROR = 0;
	}
	//--------------------------------------PID误差限幅
	if(VF32_TEM_PI_ERROR > 0.2f)
		VF32_TEM_PI_ERROR = 0.2f;
	else if(VF32_TEM_PI_ERROR < -0.2f)
		VF32_TEM_PI_ERROR = -0.2f;
	//--------------------------------------获取当前误差的增量
	VF32_TEM_PI_ERROR_Delta = VF32_TEM_PI_ERROR - *Pid->error_Last;
	//--------------------------------------获取当前误差增量的导数
	VF32_TEM_PI_pERROR_Delta = VF32_TEM_PI_ERROR_Delta - *Pid->Delta_Last;
	//--------------------------------------采用梯形积分法求取积分值
	VF32_TEM_TEMP = *Pid->KI * (VF32_TEM_PI_ERROR + *Pid->error_Last)/* - (*Pid->PI_OUT - *Pid->Regulator_out)*/;
	//--------------------------------------PID运算
	*Pid->PI_OUT += (*Pid->KP * VF32_TEM_PI_ERROR_Delta + VF32_TEM_TEMP + *Pid->KD * VF32_TEM_PI_pERROR_Delta + *Pid->Delta_Cur * *Pid->KC);
	//--------------------------------------
	if(*Pid->PI_OUT > Pid->UpperLimit)
		*Pid->PI_OUT = Pid->UpperLimit;
	else if(*Pid->PI_OUT < Pid->LowerLimit)
		*Pid->PI_OUT = Pid->LowerLimit;
	//--------------------------------------存储当前的误差和上次的误差
	*Pid->Delta_Last = VF32_TEM_PI_ERROR_Delta;
	*Pid->error_Last = VF32_TEM_PI_ERROR;
	//--------------------------------------PID输出
	return *Pid->PI_OUT;
}
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
static real32 PIDICF_INV1(pidifv1 *Pid)
{
	real32	VF32_TEM_PI_ERROR;				/*PID运算的当前误差值*/
	real32	VF32_TEM_PI_ERROR_Delta;		/*PID运算的误差值增量*/
	real32	VF32_TEM_PI_pERROR_Delta;		/*PID运算的误差值增量的导数*/
	real32	VF32_TEM_TEMP;
	//--------------------------------------获取当前误差
	VF32_TEM_TEMP = (*Pid->current_Value - *Pid->target_Value);
	//--------------------------------------禁带计算
	if(*Pid->target_Value < *Pid->KFuzz * 1.2f)//当给定值小于1.2f倍禁带宽度时，不对禁带宽度进行干预
	{
		VF32_TEM_PI_ERROR = VF32_TEM_TEMP * *Pid->ReferenceRatio;
	}
	else
	{
		if(VF32_TEM_TEMP >= *Pid->KFuzz)
			VF32_TEM_PI_ERROR = (VF32_TEM_TEMP - *Pid->KFuzz) * *Pid->ReferenceRatio;
		else if(VF32_TEM_TEMP <= -(*Pid->KFuzz))
			VF32_TEM_PI_ERROR = (VF32_TEM_TEMP + *Pid->KFuzz) * *Pid->ReferenceRatio;
		else
			VF32_TEM_PI_ERROR = 0;
	}
	//--------------------------------------PID误差限幅
	if(VF32_TEM_PI_ERROR > 0.02f)
		VF32_TEM_PI_ERROR = 0.02f;
	else if(VF32_TEM_PI_ERROR < -0.02f)
		VF32_TEM_PI_ERROR = -0.02f;
	//--------------------------------------获取当前误差的增量
	VF32_TEM_PI_ERROR_Delta = VF32_TEM_PI_ERROR - *Pid->error_Last;
	//--------------------------------------获取当前误差增量的导数
	VF32_TEM_PI_pERROR_Delta = VF32_TEM_PI_ERROR_Delta - *Pid->Delta_Last;
	//--------------------------------------采用梯形积分法求取积分值
	VF32_TEM_TEMP = *Pid->KI * (VF32_TEM_PI_ERROR + *Pid->error_Last)/* - (*Pid->PI_OUT - *Pid->Regulator_out)*/;
	//--------------------------------------PID运算
	*Pid->PI_OUT += (*Pid->KP * VF32_TEM_PI_ERROR_Delta + VF32_TEM_TEMP + *Pid->KD * VF32_TEM_PI_pERROR_Delta + *Pid->Delta_Cur * *Pid->KC);
	//--------------------------------------
	if(*Pid->PI_OUT > Pid->UpperLimit)
		*Pid->PI_OUT = Pid->UpperLimit;
	else if(*Pid->PI_OUT < Pid->LowerLimit)
		*Pid->PI_OUT = Pid->LowerLimit;
	//--------------------------------------存储当前的误差和上次的误差
	*Pid->Delta_Last = VF32_TEM_PI_ERROR_Delta;
	*Pid->error_Last = VF32_TEM_PI_ERROR;
	//--------------------------------------PID输出
	return *Pid->PI_OUT;
}
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
static real32 PIDICFV2(pidifv2 *Pid)
{
	real32	VF32_TEM_PI_ERROR;				/*PID运算的当前误差值*/
	real32	VF32_TEM_PI_ERROR_Delta;		/*PID运算的误差值增量*/
	real32	VF32_TEM_PI_pERROR_Delta;		/*PID运算的误差值增量的导数*/
	real32	VF32_TEM_TEMP;
	//--------------------------------------获取当前误差
	VF32_TEM_TEMP = (*Pid->target_Value - *Pid->current_Value);
	//--------------------------------------禁带计算
	if(*Pid->target_Value < *Pid->KFuzz * 1.2f)//当给定值小于1.2f倍禁带宽度时，不对禁带宽度进行干预
	{
		VF32_TEM_PI_ERROR = VF32_TEM_TEMP * *Pid->ReferenceRatio;
	}
	else
	{
		if(VF32_TEM_TEMP >= *Pid->KFuzz)
			VF32_TEM_PI_ERROR = (VF32_TEM_TEMP - *Pid->KFuzz) * *Pid->ReferenceRatio;
		else if(VF32_TEM_TEMP <= -(*Pid->KFuzz))
			VF32_TEM_PI_ERROR = (VF32_TEM_TEMP + *Pid->KFuzz) * *Pid->ReferenceRatio;
		else
			VF32_TEM_PI_ERROR = 0;
	}
	//--------------------------------------PID误差限幅
	if(VF32_TEM_PI_ERROR > Pid->ErrorLimit)
		VF32_TEM_PI_ERROR = Pid->ErrorLimit;
	else if(VF32_TEM_PI_ERROR < -Pid->ErrorLimit)
		VF32_TEM_PI_ERROR = -Pid->ErrorLimit;
	//--------------------------------------获取当前误差的增量
	VF32_TEM_PI_ERROR_Delta = VF32_TEM_PI_ERROR - *Pid->error_Last;
	//--------------------------------------获取当前误差增量的导数
	VF32_TEM_PI_pERROR_Delta = VF32_TEM_PI_ERROR_Delta - *Pid->Delta_Last;
	//--------------------------------------采用梯形积分法求取积分值
	VF32_TEM_TEMP = *Pid->KI * (VF32_TEM_PI_ERROR + *Pid->error_Last)/* - (*Pid->PI_OUT - *Pid->Regulator_out)*/;
	//--------------------------------------PID运算
	*Pid->PI_OUT += (*Pid->KP * VF32_TEM_PI_ERROR_Delta + VF32_TEM_TEMP + *Pid->KD * VF32_TEM_PI_pERROR_Delta + *Pid->Delta_Cur * *Pid->KC);
	//--------------------------------------
	if(*Pid->PI_OUT > Pid->UpperLimit)
		*Pid->PI_OUT = Pid->UpperLimit;
	else if(*Pid->PI_OUT < Pid->LowerLimit)
		*Pid->PI_OUT = Pid->LowerLimit;
	//--------------------------------------存储当前的误差和上次的误差
	*Pid->Delta_Last = VF32_TEM_PI_ERROR_Delta;
	*Pid->error_Last = VF32_TEM_PI_ERROR;
	//--------------------------------------PID输出
	return *Pid->PI_OUT;
}
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
static real32 PIDICF_INV2(pidifv2 *Pid)
{
	real32	VF32_TEM_PI_ERROR;				/*PID运算的当前误差值*/
	real32	VF32_TEM_PI_ERROR_Delta;		/*PID运算的误差值增量*/
	real32	VF32_TEM_PI_pERROR_Delta;		/*PID运算的误差值增量的导数*/
	real32	VF32_TEM_TEMP;
	//--------------------------------------获取当前误差
	VF32_TEM_TEMP = (*Pid->current_Value - *Pid->target_Value);
	//--------------------------------------禁带计算
	if(*Pid->target_Value < *Pid->KFuzz * 1.2f)//当给定值小于1.2f倍禁带宽度时，不对禁带宽度进行干预
	{
		VF32_TEM_PI_ERROR = VF32_TEM_TEMP * *Pid->ReferenceRatio;
	}
	else
	{
		if(VF32_TEM_TEMP >= *Pid->KFuzz)
			VF32_TEM_PI_ERROR = (VF32_TEM_TEMP - *Pid->KFuzz) * *Pid->ReferenceRatio;
		else if(VF32_TEM_TEMP <= -(*Pid->KFuzz))
			VF32_TEM_PI_ERROR = (VF32_TEM_TEMP + *Pid->KFuzz) * *Pid->ReferenceRatio;
		else
			VF32_TEM_PI_ERROR = 0;
	}
	//--------------------------------------PID误差限幅
	if(VF32_TEM_PI_ERROR > Pid->ErrorLimit)
		VF32_TEM_PI_ERROR = Pid->ErrorLimit;
	else if(VF32_TEM_PI_ERROR < -Pid->ErrorLimit)
		VF32_TEM_PI_ERROR = -Pid->ErrorLimit;
	//--------------------------------------获取当前误差的增量
	VF32_TEM_PI_ERROR_Delta = VF32_TEM_PI_ERROR - *Pid->error_Last;
	//--------------------------------------获取当前误差增量的导数
	VF32_TEM_PI_pERROR_Delta = VF32_TEM_PI_ERROR_Delta - *Pid->Delta_Last;
	//--------------------------------------采用梯形积分法求取积分值
	VF32_TEM_TEMP = *Pid->KI * (VF32_TEM_PI_ERROR + *Pid->error_Last)/* - (*Pid->PI_OUT - *Pid->Regulator_out)*/;
	//--------------------------------------PID运算
	*Pid->PI_OUT += (*Pid->KP * VF32_TEM_PI_ERROR_Delta + VF32_TEM_TEMP + *Pid->KD * VF32_TEM_PI_pERROR_Delta + *Pid->Delta_Cur * *Pid->KC);
	//--------------------------------------
	if(*Pid->PI_OUT > Pid->UpperLimit)
		*Pid->PI_OUT = Pid->UpperLimit;
	else if(*Pid->PI_OUT < Pid->LowerLimit)
		*Pid->PI_OUT = Pid->LowerLimit;
	//--------------------------------------存储当前的误差和上次的误差
	*Pid->Delta_Last = VF32_TEM_PI_ERROR_Delta;
	*Pid->error_Last = VF32_TEM_PI_ERROR;
	//--------------------------------------PID输出
	return *Pid->PI_OUT;
}
/*
 * Copyright 长沙奥托自动化技术有限公司
 * All rights reserved.
 *
 * 函数名称：float my_sin(float Thita)
 * *
 * 摘    要：实现单精度浮点型弧度数据正弦运算，返回单精度浮点型的结果
 *
 * 当前版本：1.0
 * 作    者：邓王飞
 * 完成日期：2012年12月2日
 * 备注：角度范围：-∞~+∞,角度在-2π~2π之间，执行一次运算占用1.072uS，NIOS II主频为125MHZ
 */

static real32 os_sinf(real32 radian)
{
	real32 y = 0,y1 = 0;
	while(radian > 3.141593653f)
	radian -= 6.2831853f;
	while(radian < -3.141593653f)
	radian += 6.2831853f;
	y1 = radian * radian;
	y = ((((2.169938182926827e-06f * y1 - 1.931004385730262e-04f) * y1 + 0.008312007531852f) * y1 - 0.166631756022613f) * y1 + 0.999984114910519f) * radian;
	return y;
}
/*
 * Copyright 长沙奥托自动化技术有限公司
 * All rights reserved.
 *
 * 函数名称：float my_cos(float Thita)
 * *
 * 摘    要：实现单精度浮点型弧度数据余弦运算，返回单精度浮点型的结果
 *
 * 当前版本：1.0
 * 作    者：邓王飞
 * 完成日期：2012年12月2日
 * 备注：角度范围：-∞~+∞,角度在-2π~2π之间，执行一次运算占用1.277uS，NIOS II主频为125MHZ
 */

static real32 os_cosf(real32 radian)
{
	real32 y = 0,y1 = 0;
	while(radian > 3.141593653f)
	radian -= 6.2831853f;
	while(radian < -3.141593653f)
	radian += 6.2831853f;
	y1 = radian * radian;
	y = ((((-2.216328888405680e-07f * y1 + 1.903289177264872e-05f) * y1 - 0.001343579731697f) * y1 + 0.041519694923020f) * y1 - 0.499833620946754f) * y1 + 0.999970198701940f;
	return y;
}
volatile htos_arith HTOS_Arith = {
		CLARKE,
		SYNC_RMS_V1,
		os_cosf,
		os_sinf,
		PID,
		PIDI,
		PIDIC,
		PIDICF,
		PIDICF_INV,
		PIDICFV1,
		PIDICF_INV1,
		PIDICFV2,
		PIDICF_INV2
};
