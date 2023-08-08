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
 * File: Typedef.h
 * Author: Hunter.ORG
 * Created on: 2016年2月16日
 *
 * -Modified for F28069 DSP data types for int16 and int32 Joe Zhou 01-03-2018
 *
 * Revision history: 1.1
 */

/*
 * This is a guard condition so that contents of this file are not included
 * more than once.
 */


#ifndef TYPEDEF_H_
#define TYPEDEF_H_

typedef unsigned char		uint8;		/*无符号8位*/
typedef char				int8;		/*有符号8位*/
typedef unsigned short		uint16;		/*无符号16位*/
typedef short				int16;		//
typedef unsigned int		uint32;		/*无符号32位*/
typedef int					int32;		//
typedef unsigned long long	uint64;		/*无符号64位*/
typedef long long			int64;		/*又符号64位*/
typedef float				real32;		/*单精度浮点型*/
typedef	double				real64;		/*双精度浮点型*/
//----------------------------------------------特殊数据类型定义
typedef struct
{
	real32 alpha;
	real32 beta;
}complex;								/*复数*/
#endif /* TYPEDEF_H_ */
