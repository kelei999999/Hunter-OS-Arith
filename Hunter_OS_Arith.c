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
 * Created on: 2014��4��12��
 * Revision history: 1.1
 */

/*
 * This is a guard condition so that contents of this file are not included
 * more than once.
 */
#include "Hunter_OS_Arith.h"
#include <math.h>
//****************************************************************************
// @������			complex CLARKE(real32_T A,real32_T B,real32_T C)
//----------------------------------------------------------------------------
// @����			����ʵ����ckarke�ȷ��任
//
//----------------------------------------------------------------------------
// @����			���ྲֹ�����µ�˲ʱֵ
//
//----------------------------------------------------------------------------
// @���			���ྲֹ�����µ�˲ʱֵ
//
//----------------------------------------------------------------------------
// @����ֵ			complex
//					complex����Ϊ:
//						typedef struct
//						{
//							float alafa;
//							float beita;
//						}complex;
//
//----------------------------------------------------------------------------
// @����			2017��3��21��
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
// @������			real32 PID(pid *Pid)
//----------------------------------------------------------------------------
// @����			����ʵ����PID����
//					���ERROR = Ŀ��ֵ - ��ǰֵ
//----------------------------------------------------------------------------
// @����			pid���������pid�ṹ��
//
//----------------------------------------------------------------------------
// @���			PID������ֵ
//
//----------------------------------------------------------------------------
// @����ֵ			real32 PID������ֵ
//----------------------------------------------------------------------------
// @����			2017��4��20��
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
// @������			real32 PIDI(pidi *Pid)
//----------------------------------------------------------------------------
// @����			����ʵ��������ʽPID����
//					���ERROR = Ŀ��ֵ - ��ǰֵ
//----------------------------------------------------------------------------
// @����			pid���������pidi�ṹ��
//
//----------------------------------------------------------------------------
// @���			PID������ֵ
//
//----------------------------------------------------------------------------
// @����ֵ			real32 PID������ֵ
//----------------------------------------------------------------------------
// @����			2017��4��20��
//
//****************************************************************************
//extern real32 					VF32_TEM_TEMP1,VF32_TEM_TEMP2,VF32_TEM_TEMP3;
static real32 PIDI(pidi *Pid)
{
	real32	VF32_TEM_PI_ERROR;				/*PID����ĵ�ǰ���ֵ*/
	real32	VF32_TEM_PI_ERROR_Delta;		/*PID��������ֵ����*/
	real32	VF32_TEM_PI_pERROR_Delta;		/*PID��������ֵ�����ĵ���*/
	real32	PIDI_TEMP;
	//--------------------------------------��ȡ��ǰ���
	VF32_TEM_PI_ERROR = (*Pid->target_Value - *Pid->current_Value) * *Pid->ReferenceRatio;
	//--------------------------------------PID����޷�
	if(VF32_TEM_PI_ERROR > 0.2f)
		VF32_TEM_PI_ERROR = 0.2f;
	else if(VF32_TEM_PI_ERROR < -0.2f)
		VF32_TEM_PI_ERROR = -0.2f;
	//--------------------------------------��ȡ��ǰ��������
	VF32_TEM_PI_ERROR_Delta = VF32_TEM_PI_ERROR - *Pid->error_Last;
	//--------------------------------------��ȡ��ǰ��������ĵ���
	VF32_TEM_PI_pERROR_Delta = VF32_TEM_PI_ERROR_Delta - *Pid->Delta_Last;
	//--------------------------------------�������λ��ַ���ȡ����ֵ
	PIDI_TEMP = *Pid->KI * (VF32_TEM_PI_ERROR + *Pid->error_Last) - (*Pid->PI_OUT - *Pid->Regulator_out);
	//--------------------------------------PID����
	*Pid->PI_OUT += (*Pid->KP * VF32_TEM_PI_ERROR_Delta + PIDI_TEMP + *Pid->KD * VF32_TEM_PI_pERROR_Delta);
	//--------------------------------------�洢��ǰ�������ϴε����
	*Pid->Delta_Last = VF32_TEM_PI_ERROR_Delta;
	*Pid->error_Last = VF32_TEM_PI_ERROR;
	//--------------------------------------PID���
	return *Pid->PI_OUT;
}
//****************************************************************************
// @������			real32 PIDC(pidc *Pid)
//----------------------------------------------------------------------------
// @����			����ʵ���˴�������ֲ���������ʽPID����
//					���ERROR = Ŀ��ֵ - ��ǰֵ
//----------------------------------------------------------------------------
// @����			pid���������pidc�ṹ��
//
//----------------------------------------------------------------------------
// @���			PID������ֵ
//
//----------------------------------------------------------------------------
// @����ֵ			real32 PID������ֵ
//----------------------------------------------------------------------------
// @����			2017��4��20��
//
//****************************************************************************
//extern real32 					VF32_TEM_TEMP1,VF32_TEM_TEMP2,VF32_TEM_TEMP3;
static real32 PIDIC(pidic *Pid)
{
	real32	VF32_TEM_PI_ERROR;				/*PID����ĵ�ǰ���ֵ*/
	real32	VF32_TEM_PI_ERROR_Delta;		/*PID��������ֵ����*/
	real32	VF32_TEM_PI_pERROR_Delta;		/*PID��������ֵ�����ĵ���*/
	real32	PIDI_TEMP;
	//--------------------------------------��ȡ��ǰ���
	VF32_TEM_PI_ERROR = (*Pid->target_Value - *Pid->current_Value) * *Pid->ReferenceRatio;
	//--------------------------------------PID����޷�
	if(VF32_TEM_PI_ERROR > 0.8f)
		VF32_TEM_PI_ERROR = 0.8f;
	else if(VF32_TEM_PI_ERROR < -0.8f)
		VF32_TEM_PI_ERROR = -0.8f;
	//--------------------------------------��ȡ��ǰ��������
	VF32_TEM_PI_ERROR_Delta = VF32_TEM_PI_ERROR - *Pid->error_Last;
	//--------------------------------------��ȡ��ǰ��������ĵ���
	VF32_TEM_PI_pERROR_Delta = VF32_TEM_PI_ERROR_Delta - *Pid->Delta_Last;
	//--------------------------------------�������λ��ַ���ȡ����ֵ
	PIDI_TEMP = *Pid->KI * (VF32_TEM_PI_ERROR + *Pid->error_Last) - (*Pid->PI_OUT - *Pid->Regulator_out);
	//--------------------------------------PID����
	*Pid->PI_OUT += (*Pid->KP * VF32_TEM_PI_ERROR_Delta + PIDI_TEMP + *Pid->KD * VF32_TEM_PI_pERROR_Delta + *Pid->Delta_Cur * *Pid->KC);
	//--------------------------------------�洢��ǰ�������ϴε����
	*Pid->Delta_Last = VF32_TEM_PI_ERROR_Delta;
	*Pid->error_Last = VF32_TEM_PI_ERROR;
	//--------------------------------------PID���
	return *Pid->PI_OUT;
}
//****************************************************************************
// @������			real32 PIDICF(pidif *Pid)
//----------------------------------------------------------------------------
// @����			����ʵ���˴�������ֲ�����ģ��������ȵ�����ʽPID����
//					������ֵС��2���������ʱ��������Ƚ���������
//					���ERROR = Ŀ��ֵ - ��ǰֵ
//----------------------------------------------------------------------------
// @����			pid���������pidif�ṹ��
//
//----------------------------------------------------------------------------
// @���			PID������ֵ
//
//----------------------------------------------------------------------------
// @����ֵ			real32 PID������ֵ
//----------------------------------------------------------------------------
// @����			2017��4��20��
//
//****************************************************************************
static real32 PIDICF(pidif *Pid)
{
	real32	VF32_TEM_PI_ERROR;				/*PID����ĵ�ǰ���ֵ*/
	real32	VF32_TEM_PI_ERROR_Delta;		/*PID��������ֵ����*/
	real32	VF32_TEM_PI_pERROR_Delta;		/*PID��������ֵ�����ĵ���*/
	real32	VF32_TEM_TEMP;
	//--------------------------------------��ȡ��ǰ���
	VF32_TEM_TEMP = (*Pid->target_Value - *Pid->current_Value);
	//--------------------------------------��������
	if(*Pid->target_Value < *Pid->KFuzz * 1.2f)//������ֵС��1.2f���������ʱ�����Խ�����Ƚ��и�Ԥ
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
	//--------------------------------------PID����޷�
	if(VF32_TEM_PI_ERROR > 0.4f)
		VF32_TEM_PI_ERROR = 0.4f;
	else if(VF32_TEM_PI_ERROR < -0.4f)
		VF32_TEM_PI_ERROR = -0.4f;
	//--------------------------------------��ȡ��ǰ��������
	VF32_TEM_PI_ERROR_Delta = VF32_TEM_PI_ERROR - *Pid->error_Last;
	//--------------------------------------��ȡ��ǰ��������ĵ���
	VF32_TEM_PI_pERROR_Delta = VF32_TEM_PI_ERROR_Delta - *Pid->Delta_Last;
	//--------------------------------------�������λ��ַ���ȡ����ֵ
	VF32_TEM_TEMP = *Pid->KI * (VF32_TEM_PI_ERROR + *Pid->error_Last) - (*Pid->PI_OUT - *Pid->Regulator_out);
	//--------------------------------------PID����
	*Pid->PI_OUT += (*Pid->KP * VF32_TEM_PI_ERROR_Delta + VF32_TEM_TEMP + *Pid->KD * VF32_TEM_PI_pERROR_Delta + *Pid->Delta_Cur * *Pid->KC);
	//--------------------------------------�洢��ǰ�������ϴε����
	*Pid->Delta_Last = VF32_TEM_PI_ERROR_Delta;
	*Pid->error_Last = VF32_TEM_PI_ERROR;
	//--------------------------------------PID���
	return *Pid->PI_OUT;
}
//****************************************************************************
// @������			real32 PIDICF_INV(pidif *Pid)
//----------------------------------------------------------------------------
// @����			����ʵ���˴�������ֲ�����ģ��������ȵ�����ʽ����PID����
//					������ֵС��2���������ʱ��������Ƚ���������
//					���ERROR = ��ǰֵ - Ŀ��ֵ
//----------------------------------------------------------------------------
// @����			pid���������pidif�ṹ��
//
//----------------------------------------------------------------------------
// @���			PID������ֵ
//
//----------------------------------------------------------------------------
// @����ֵ			real32 PID������ֵ
//----------------------------------------------------------------------------
// @����			2017��4��26��
//
//****************************************************************************
static real32 PIDICF_INV(pidif *Pid)
{
	real32	VF32_TEM_PI_ERROR;				/*PID����ĵ�ǰ���ֵ*/
	real32	VF32_TEM_PI_ERROR_Delta;		/*PID��������ֵ����*/
	real32	VF32_TEM_PI_pERROR_Delta;		/*PID��������ֵ�����ĵ���*/
	real32	VF32_TEM_TEMP;
	//--------------------------------------��ȡ��ǰ���
	VF32_TEM_TEMP = (*Pid->current_Value - *Pid->target_Value);
	//--------------------------------------��������
	if(*Pid->target_Value < *Pid->KFuzz * 1.2f)//������ֵС��1.2f���������ʱ�����Խ�����Ƚ��и�Ԥ
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
	//--------------------------------------PID����޷�
	if(VF32_TEM_PI_ERROR > 0.4f)
		VF32_TEM_PI_ERROR = 0.4f;
	else if(VF32_TEM_PI_ERROR < -0.4f)
		VF32_TEM_PI_ERROR = -0.4f;
	//--------------------------------------��ȡ��ǰ��������
	VF32_TEM_PI_ERROR_Delta = VF32_TEM_PI_ERROR - *Pid->error_Last;
	//--------------------------------------��ȡ��ǰ��������ĵ���
	VF32_TEM_PI_pERROR_Delta = VF32_TEM_PI_ERROR_Delta - *Pid->Delta_Last;
	//--------------------------------------�������λ��ַ���ȡ����ֵ
	VF32_TEM_TEMP = *Pid->KI * (VF32_TEM_PI_ERROR + *Pid->error_Last) - (*Pid->PI_OUT - *Pid->Regulator_out);
	//--------------------------------------PID����
	*Pid->PI_OUT += (*Pid->KP * VF32_TEM_PI_ERROR_Delta + VF32_TEM_TEMP + *Pid->KD * VF32_TEM_PI_pERROR_Delta + *Pid->Delta_Cur * *Pid->KC);
	//--------------------------------------�洢��ǰ�������ϴε����
	*Pid->Delta_Last = VF32_TEM_PI_ERROR_Delta;
	*Pid->error_Last = VF32_TEM_PI_ERROR;
	//--------------------------------------PID���
	return *Pid->PI_OUT;
}
//****************************************************************************
// @������		real32 SYNC_RMS_V1(rms_process *RMS)
//----------------------------------------------------------------------------
// @����                       ����ʵ����αͬ��RMS����
//
//----------------------------------------------------------------------------
// @����			rms_process *RMS:���rms_process�ṹ��
//
//----------------------------------------------------------------------------
// @���
//
//----------------------------------------------------------------------------
// @����ֵ		RMS��ֵ
//----------------------------------------------------------------------------
// @����                       2017��5��28��
//
//****************************************************************************
static real32 SYNC_RMS_V1(rms_process *RMS)
{
	real32 temp = *RMS->input.INS;
	//--------------------------------------------�ҹ����
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
	//-------------------------------------------����ƽ���ۺ�
	RMS->temp.RMS_ACC += (temp * temp);
	//-------------------------------------------��¼�ۺʹ���
	RMS->temp.TIME += 1;
	//-------------------------------------------
	if(RMS->temp.TIME >= RMS->input.Count_MAX)
		RMS->temp.ZERO |= 1;
	//-------------------------------------------���������
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
// @������			real32 PIDICFV1(pidifv1 *Pid)
//----------------------------------------------------------------------------
// @����			����ʵ���˴�������ֲ�����ģ��������ȵ�����ʽPID����
//					������ֵС��2���������ʱ��������Ƚ���������
//					���ERROR = Ŀ��ֵ - ��ǰֵ
//----------------------------------------------------------------------------
// @����			pid���������pidifv1�ṹ��
//
//----------------------------------------------------------------------------
// @���			PID������ֵ
//
//----------------------------------------------------------------------------
// @����ֵ			real32 PID������ֵ
//----------------------------------------------------------------------------
// @����			2017��4��20��
//
//****************************************************************************
static real32 PIDICFV1(pidifv1 *Pid)
{
	real32	VF32_TEM_PI_ERROR;				/*PID����ĵ�ǰ���ֵ*/
	real32	VF32_TEM_PI_ERROR_Delta;		/*PID��������ֵ����*/
	real32	VF32_TEM_PI_pERROR_Delta;		/*PID��������ֵ�����ĵ���*/
	real32	VF32_TEM_TEMP;
	//--------------------------------------��ȡ��ǰ���
	VF32_TEM_TEMP = (*Pid->target_Value - *Pid->current_Value);
	//--------------------------------------��������
	if(*Pid->target_Value < *Pid->KFuzz * 1.2f)//������ֵС��1.2f���������ʱ�����Խ�����Ƚ��и�Ԥ
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
	//--------------------------------------PID����޷�
	if(VF32_TEM_PI_ERROR > 0.2f)
		VF32_TEM_PI_ERROR = 0.2f;
	else if(VF32_TEM_PI_ERROR < -0.2f)
		VF32_TEM_PI_ERROR = -0.2f;
	//--------------------------------------��ȡ��ǰ��������
	VF32_TEM_PI_ERROR_Delta = VF32_TEM_PI_ERROR - *Pid->error_Last;
	//--------------------------------------��ȡ��ǰ��������ĵ���
	VF32_TEM_PI_pERROR_Delta = VF32_TEM_PI_ERROR_Delta - *Pid->Delta_Last;
	//--------------------------------------�������λ��ַ���ȡ����ֵ
	VF32_TEM_TEMP = *Pid->KI * (VF32_TEM_PI_ERROR + *Pid->error_Last)/* - (*Pid->PI_OUT - *Pid->Regulator_out)*/;
	//--------------------------------------PID����
	*Pid->PI_OUT += (*Pid->KP * VF32_TEM_PI_ERROR_Delta + VF32_TEM_TEMP + *Pid->KD * VF32_TEM_PI_pERROR_Delta + *Pid->Delta_Cur * *Pid->KC);
	//--------------------------------------
	if(*Pid->PI_OUT > Pid->UpperLimit)
		*Pid->PI_OUT = Pid->UpperLimit;
	else if(*Pid->PI_OUT < Pid->LowerLimit)
		*Pid->PI_OUT = Pid->LowerLimit;
	//--------------------------------------�洢��ǰ�������ϴε����
	*Pid->Delta_Last = VF32_TEM_PI_ERROR_Delta;
	*Pid->error_Last = VF32_TEM_PI_ERROR;
	//--------------------------------------PID���
	return *Pid->PI_OUT;
}
//****************************************************************************
// @������			real32 PIDICF_INV1(pidifv1 *Pid)
//----------------------------------------------------------------------------
// @����			����ʵ���˴�������ֲ�����ģ��������ȵ�����ʽ����PID����
//					������ֵС��2���������ʱ��������Ƚ���������
//					���ERROR = ��ǰֵ - Ŀ��ֵ
//----------------------------------------------------------------------------
// @����			pid���������pidifv1�ṹ��
//
//----------------------------------------------------------------------------
// @���			PID������ֵ
//
//----------------------------------------------------------------------------
// @����ֵ			real32 PID������ֵ
//----------------------------------------------------------------------------
// @����			2017��4��26��
//
//****************************************************************************
static real32 PIDICF_INV1(pidifv1 *Pid)
{
	real32	VF32_TEM_PI_ERROR;				/*PID����ĵ�ǰ���ֵ*/
	real32	VF32_TEM_PI_ERROR_Delta;		/*PID��������ֵ����*/
	real32	VF32_TEM_PI_pERROR_Delta;		/*PID��������ֵ�����ĵ���*/
	real32	VF32_TEM_TEMP;
	//--------------------------------------��ȡ��ǰ���
	VF32_TEM_TEMP = (*Pid->current_Value - *Pid->target_Value);
	//--------------------------------------��������
	if(*Pid->target_Value < *Pid->KFuzz * 1.2f)//������ֵС��1.2f���������ʱ�����Խ�����Ƚ��и�Ԥ
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
	//--------------------------------------PID����޷�
	if(VF32_TEM_PI_ERROR > 0.02f)
		VF32_TEM_PI_ERROR = 0.02f;
	else if(VF32_TEM_PI_ERROR < -0.02f)
		VF32_TEM_PI_ERROR = -0.02f;
	//--------------------------------------��ȡ��ǰ��������
	VF32_TEM_PI_ERROR_Delta = VF32_TEM_PI_ERROR - *Pid->error_Last;
	//--------------------------------------��ȡ��ǰ��������ĵ���
	VF32_TEM_PI_pERROR_Delta = VF32_TEM_PI_ERROR_Delta - *Pid->Delta_Last;
	//--------------------------------------�������λ��ַ���ȡ����ֵ
	VF32_TEM_TEMP = *Pid->KI * (VF32_TEM_PI_ERROR + *Pid->error_Last)/* - (*Pid->PI_OUT - *Pid->Regulator_out)*/;
	//--------------------------------------PID����
	*Pid->PI_OUT += (*Pid->KP * VF32_TEM_PI_ERROR_Delta + VF32_TEM_TEMP + *Pid->KD * VF32_TEM_PI_pERROR_Delta + *Pid->Delta_Cur * *Pid->KC);
	//--------------------------------------
	if(*Pid->PI_OUT > Pid->UpperLimit)
		*Pid->PI_OUT = Pid->UpperLimit;
	else if(*Pid->PI_OUT < Pid->LowerLimit)
		*Pid->PI_OUT = Pid->LowerLimit;
	//--------------------------------------�洢��ǰ�������ϴε����
	*Pid->Delta_Last = VF32_TEM_PI_ERROR_Delta;
	*Pid->error_Last = VF32_TEM_PI_ERROR;
	//--------------------------------------PID���
	return *Pid->PI_OUT;
}
//****************************************************************************
// @������			real32 PIDICFV2(pidifv2 *Pid)
//----------------------------------------------------------------------------
// @����			����ʵ���˴�������ֲ�����ģ��������ȵ�����ʽPID����
//					������ֵС��2���������ʱ��������Ƚ���������
//					���ERROR = Ŀ��ֵ - ��ǰֵ
//----------------------------------------------------------------------------
// @����			pid���������pidifv2�ṹ��
//
//----------------------------------------------------------------------------
// @���			PID������ֵ
//
//----------------------------------------------------------------------------
// @����ֵ			real32 PID������ֵ
//----------------------------------------------------------------------------
// @����			2017��4��20��
//
//****************************************************************************
static real32 PIDICFV2(pidifv2 *Pid)
{
	real32	VF32_TEM_PI_ERROR;				/*PID����ĵ�ǰ���ֵ*/
	real32	VF32_TEM_PI_ERROR_Delta;		/*PID��������ֵ����*/
	real32	VF32_TEM_PI_pERROR_Delta;		/*PID��������ֵ�����ĵ���*/
	real32	VF32_TEM_TEMP;
	//--------------------------------------��ȡ��ǰ���
	VF32_TEM_TEMP = (*Pid->target_Value - *Pid->current_Value);
	//--------------------------------------��������
	if(*Pid->target_Value < *Pid->KFuzz * 1.2f)//������ֵС��1.2f���������ʱ�����Խ�����Ƚ��и�Ԥ
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
	//--------------------------------------PID����޷�
	if(VF32_TEM_PI_ERROR > Pid->ErrorLimit)
		VF32_TEM_PI_ERROR = Pid->ErrorLimit;
	else if(VF32_TEM_PI_ERROR < -Pid->ErrorLimit)
		VF32_TEM_PI_ERROR = -Pid->ErrorLimit;
	//--------------------------------------��ȡ��ǰ��������
	VF32_TEM_PI_ERROR_Delta = VF32_TEM_PI_ERROR - *Pid->error_Last;
	//--------------------------------------��ȡ��ǰ��������ĵ���
	VF32_TEM_PI_pERROR_Delta = VF32_TEM_PI_ERROR_Delta - *Pid->Delta_Last;
	//--------------------------------------�������λ��ַ���ȡ����ֵ
	VF32_TEM_TEMP = *Pid->KI * (VF32_TEM_PI_ERROR + *Pid->error_Last)/* - (*Pid->PI_OUT - *Pid->Regulator_out)*/;
	//--------------------------------------PID����
	*Pid->PI_OUT += (*Pid->KP * VF32_TEM_PI_ERROR_Delta + VF32_TEM_TEMP + *Pid->KD * VF32_TEM_PI_pERROR_Delta + *Pid->Delta_Cur * *Pid->KC);
	//--------------------------------------
	if(*Pid->PI_OUT > Pid->UpperLimit)
		*Pid->PI_OUT = Pid->UpperLimit;
	else if(*Pid->PI_OUT < Pid->LowerLimit)
		*Pid->PI_OUT = Pid->LowerLimit;
	//--------------------------------------�洢��ǰ�������ϴε����
	*Pid->Delta_Last = VF32_TEM_PI_ERROR_Delta;
	*Pid->error_Last = VF32_TEM_PI_ERROR;
	//--------------------------------------PID���
	return *Pid->PI_OUT;
}
//****************************************************************************
// @������			real32 PIDICF_INV2(pidifv2 *Pid)
//----------------------------------------------------------------------------
// @����			����ʵ���˴�������ֲ�����ģ��������ȵ�����ʽ����PID����
//					������ֵС��2���������ʱ��������Ƚ���������
//					���ERROR = ��ǰֵ - Ŀ��ֵ
//----------------------------------------------------------------------------
// @����			pid���������pidifv2�ṹ��
//
//----------------------------------------------------------------------------
// @���			PID������ֵ
//
//----------------------------------------------------------------------------
// @����ֵ			real32 PID������ֵ
//----------------------------------------------------------------------------
// @����			2017��4��26��
//
//****************************************************************************
static real32 PIDICF_INV2(pidifv2 *Pid)
{
	real32	VF32_TEM_PI_ERROR;				/*PID����ĵ�ǰ���ֵ*/
	real32	VF32_TEM_PI_ERROR_Delta;		/*PID��������ֵ����*/
	real32	VF32_TEM_PI_pERROR_Delta;		/*PID��������ֵ�����ĵ���*/
	real32	VF32_TEM_TEMP;
	//--------------------------------------��ȡ��ǰ���
	VF32_TEM_TEMP = (*Pid->current_Value - *Pid->target_Value);
	//--------------------------------------��������
	if(*Pid->target_Value < *Pid->KFuzz * 1.2f)//������ֵС��1.2f���������ʱ�����Խ�����Ƚ��и�Ԥ
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
	//--------------------------------------PID����޷�
	if(VF32_TEM_PI_ERROR > Pid->ErrorLimit)
		VF32_TEM_PI_ERROR = Pid->ErrorLimit;
	else if(VF32_TEM_PI_ERROR < -Pid->ErrorLimit)
		VF32_TEM_PI_ERROR = -Pid->ErrorLimit;
	//--------------------------------------��ȡ��ǰ��������
	VF32_TEM_PI_ERROR_Delta = VF32_TEM_PI_ERROR - *Pid->error_Last;
	//--------------------------------------��ȡ��ǰ��������ĵ���
	VF32_TEM_PI_pERROR_Delta = VF32_TEM_PI_ERROR_Delta - *Pid->Delta_Last;
	//--------------------------------------�������λ��ַ���ȡ����ֵ
	VF32_TEM_TEMP = *Pid->KI * (VF32_TEM_PI_ERROR + *Pid->error_Last)/* - (*Pid->PI_OUT - *Pid->Regulator_out)*/;
	//--------------------------------------PID����
	*Pid->PI_OUT += (*Pid->KP * VF32_TEM_PI_ERROR_Delta + VF32_TEM_TEMP + *Pid->KD * VF32_TEM_PI_pERROR_Delta + *Pid->Delta_Cur * *Pid->KC);
	//--------------------------------------
	if(*Pid->PI_OUT > Pid->UpperLimit)
		*Pid->PI_OUT = Pid->UpperLimit;
	else if(*Pid->PI_OUT < Pid->LowerLimit)
		*Pid->PI_OUT = Pid->LowerLimit;
	//--------------------------------------�洢��ǰ�������ϴε����
	*Pid->Delta_Last = VF32_TEM_PI_ERROR_Delta;
	*Pid->error_Last = VF32_TEM_PI_ERROR;
	//--------------------------------------PID���
	return *Pid->PI_OUT;
}
/*
 * Copyright ��ɳ�����Զ����������޹�˾
 * All rights reserved.
 *
 * �������ƣ�float my_sin(float Thita)
 * *
 * ժ    Ҫ��ʵ�ֵ����ȸ����ͻ��������������㣬���ص����ȸ����͵Ľ��
 *
 * ��ǰ�汾��1.0
 * ��    �ߣ�������
 * ������ڣ�2012��12��2��
 * ��ע���Ƕȷ�Χ��-��~+��,�Ƕ���-2��~2��֮�䣬ִ��һ������ռ��1.072uS��NIOS II��ƵΪ125MHZ
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
 * Copyright ��ɳ�����Զ����������޹�˾
 * All rights reserved.
 *
 * �������ƣ�float my_cos(float Thita)
 * *
 * ժ    Ҫ��ʵ�ֵ����ȸ����ͻ��������������㣬���ص����ȸ����͵Ľ��
 *
 * ��ǰ�汾��1.0
 * ��    �ߣ�������
 * ������ڣ�2012��12��2��
 * ��ע���Ƕȷ�Χ��-��~+��,�Ƕ���-2��~2��֮�䣬ִ��һ������ռ��1.277uS��NIOS II��ƵΪ125MHZ
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
