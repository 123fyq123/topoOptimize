#include "singalFather.h"

SingalFather *SingalFather::instance = nullptr;

SingalFather::SingalFather() : QObject(nullptr) {}


/**************************************************
@brief   : 构造一个新对象 使用单例模式  说实话这个对象没啥用，就是为了承接 更新显示信号的东西
@author  : lee
@input   ：none
@output  ：none
@time    : none
**************************************************/
SingalFather *SingalFather::getInstance() {
	if (!instance)
		instance = new SingalFather;
	return instance;
}
