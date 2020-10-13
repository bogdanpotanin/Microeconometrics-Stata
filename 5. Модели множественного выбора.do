* Потанин Богдан Станиславович
* Микроэконометрика
* Семинар №5
* Модели множественного выбора
********
* РАЗДЕЛ №0. Экспорт данных из РМЭЗ
********
sysuse r25i_os26с.dta
* Создадим необходимые для анализа переменные
* Факт трудоустройства
generate work = 1 if uj1 <= 4
replace work = 0 if uj1 == 5
replace work = . if uj1 >= 6
* Возраст
generate age = u_age
* Посмотрим описательные статистики переменной возраста
sum age
	* Пропуски в РМЭЗ получают значения 99999997, 99999998 и 99999999,
	* поэтому, при создании переменных от них следует избавляться.
replace age = . if age > 1000
* Образование
	* Законченное высшее
generate educ = 3 if u_diplom == 6
	* Закоченное среднее специальное
replace educ = 2 if u_diplom == 5
	* Закоченное среднее
replace educ = 1 if u_diplom == 4
	* Другое
replace educ = 0 if u_diplom < 4
replace educ = . if u_diplom > 6
	* Половая принадлежность
generate male = 1 if uh5 == 1
replace male = 0 if uh5 == 2
replace male = . if uh5 > 2
	* Количество детей
generate children = uj72_172
	* Для тех индивидов, у которых нет детей, переменная на количество детей uj72_172 принимает
	* значение, равное пропуску. Из-за этого при регрессионном анализе в случае использования данной
	* переменной из выборки будут ошибочно удалены все бездетные индивиды. Во избежание этого необходимо
	* заменить пропуски на 0 для тех индивидов, которые ответили на вопрос о наличии детей отрицательно.
replace children = 0 if (uj72_171 != .) & (children == .)
	* Создадим переменные на рост, вес и BMI
gen weight = um1
replace weight = . if um1 > 10000
gen height = um2
replace height = . if um2 > 10000
gen BMI = weight / (height / 100) ^ 2
	* Создадим переменную на брак
gen marriage = 1 if u_marst == 2
replace marriage = 0 if u_marst != 2
replace marriage = . if u_marst > 100
* Посмотрим на описательные статистики наших переменных
sum work age educ male children BMI marriage
	* Отдельно для мужчин
sum work age educ male children BMI marriage if male == 1
* Можно, также, посмотреть на гистограммы и ядерные оценки плотности
histogram age
kdensity age
	* Отдельно для женщин
histogram age if male == 0
kdensity age if male == 0
********
* РАЗДЕЛ №1. Анализ симулированных данных
********
mata: mata clear
* Приступим к симуляциям
	* Количество наблюдений
mata n = 10000
	* Матрица независимых переменных
mata X = (rnormal(n, 1, 1, 2), rnormal(n, 1, 1.5, 3), rnormal(n, 1, 0, 1))
	* Случайная ошибка из распределения Гумбеля (не логистического!)
mata epsilon = (-ln(-ln(runiform(n, 1))), -ln(-ln(runiform(n, 1))), -ln(-ln(runiform(n, 1))))
	* Коэффициенты, где beta_1 - стандартизируется с целью обеспечения
	* идентифицируемости параметров модели
mata beta_1 = (0, 0, 0, 0)
mata beta_2 = (0.1, 0.2, 0.3, -0.1)
mata beta_3 = (0.2, 0.3, -0.2, 0.5)
	* Создаем латентную переменную
mata z_star_1 = beta_1[1] :+ (beta_1[2] :* X[,1]) + (beta_1[3] :* X[,2]) + (beta_1[4] :* X[,3]) :+ epsilon[,1]
mata z_star_2 = beta_2[1] :+ (beta_2[2] :* X[,1]) + (beta_2[3] :* X[,2]) + (beta_2[4] :* X[,3]) :+ epsilon[,2]
mata z_star_3 = beta_3[1] :+ (beta_3[2] :* X[,1]) + (beta_3[3] :* X[,2]) + (beta_3[4] :* X[,3]) :+ epsilon[,3]
	* Создадим наблюдаемые значения зависимой переменной
mata z = J(n,1, .)
. mata 
for(i = 1; i <= n; i++)
{
	if((z_star_2[i] <= z_star_1[i]) & (z_star_3[i] <= z_star_1[i]))
	{
		z[i,1] = 1
	}
	if((z_star_1[i] <= z_star_2[i]) & (z_star_3[i] <= z_star_2[i]))
	{
		z[i,1] = 2
	}
	if((z_star_1[i] <= z_star_3[i]) & (z_star_2[i] <= z_star_3[i]))
	{
		z[i,1] = 3
	}
}
end
* Напишем функцию правдоподобия в соответствии с представленным выше
* процессом генерации данных
. mata
void lnL(todo, x, X, z, y, g, H)
{
	beta_2_m = (x[1], x[2], x[3], x[4])
	beta_3_m = (x[5], x[6], x[7], x[8])
	
	n_m = rows(X)
	
	prob_1_nominator = exp((0))
	prob_2_nominator = exp(beta_2_m[1] :+ (beta_2_m[2] :* X[,1]) + (beta_2_m[3] :* X[,2]) + (beta_2_m[4] :* X[,3]))
	prob_3_nominator = exp(beta_3_m[1] :+ (beta_3_m[2] :* X[,1]) + (beta_3_m[3] :* X[,2]) + (beta_3_m[4] :* X[,3]))
	
	prob_denominator = prob_1_nominator :+ prob_2_nominator :+ prob_3_nominator
	
	prob_1 = prob_1_nominator :/ prob_denominator
	prob_2 = prob_2_nominator :/ prob_denominator
	prob_3 = prob_3_nominator :/ prob_denominator
	
	lnL_vector = J(n_m, 1, .)
	
	for(i = 1; i <= n_m; i++)
	{
		if(z[i,1] == 1)
		{
			lnL_vector[i, 1] = prob_1[i, 1]
		}
		if(z[i,1] == 2)
		{
			lnL_vector[i, 1] = prob_2[i, 1]
		}
		if(z[i,1] == 3)
		{
			lnL_vector[i, 1] =  prob_3[i, 1]
		}
	}
	
	y = sum(log(lnL_vector))
}
end
* Запустим оптимизатор
. mata
S = optimize_init()
optimize_init_evaluator(S, &lnL())
optimize_init_params(S, (0.15, 0.25, 0.35, -0.15, 0.25, 0.35, -0.25, 0.55))
optimize_init_argument(S, 1, X)
optimize_init_argument(S, 2, z)
optimize_init_conv_maxiter(S, 1000)
x = optimize(S)
H = optimize_result_Hessian(S)
end
* Посмотрим на полученные значения и сравним их с истинными
mata ((beta_2, beta_3) \ x)
* ЗАДАНИЯ
* 1. Проверьте, что случится, если вы неверное специфицируете количество категорий, например,
*    объединив первую и вторую в одну
* 2. Замените распределение случайной ошибки с нормального на логистическое и посмотрите, как
*    изменятся результаты оценивания параметров модели.
********
* РАЗДЕЛ №2. Приложение модели множественного выбора к реальным данным
********
* Создадим переменную на предпочитаемый тип сигарет
	* Сигареты с фильтром
gen ciga = 1 if (um74 == 2)
	* Сигареты без фильтра
replace ciga = 2 if (um74 != 2) & (um74 != .)
	* Не курит
replace ciga = 3 if (um71 == 2)
	* Посмотрим на распределение предпочтений
tab ciga
	* Построим модель для совершеннолетних мужчин
mlogit ciga c.age c.age#c.age i.educ i.marriage c.children c.BMI c.BMI#c.BMI work if(male == 1 & age >= 18)
estimates store mlogit_model_F
 * Посчитаем предельные эффекты для каждой категории
	* Средний предельный эффект
margins, dydx(*)
	* Предельные эффекты для среднего индивида
margins, dydx(*) atmeans
	* Предельные эффекты для конкретного индивида
margins, dydx(*) atmeans at(age = 50 educ = 3 BMI = 23 children = 2 marriage = 1 work = 1 wage = 50000)
* Оценим точность модели
* Получим предсказанные вероятности
predict p1 p2 p3, pr
* Получим бинарные переменные на принадлежность к каждому из классов
gen z_1 = 0 if e(sample)
gen z_2 = 0 if e(sample)
gen z_3 = 0 if e(sample)
replace z_1 = 1 if ((p1 >= p1) & (p1 >= p2) & (p1 >= p3) & e(sample))
replace z_2 = 1 if ((p2 >= p1) & (p2 >= p2) & (p2 >= p3)& e(sample))
replace z_3 = 1 if ((p3 >= p1) & (p3 >= p2) & (p3 >= p3)& e(sample))
* Создадим переменную на предсказанный класс
gen z_class = z_1 * 1 + z_2 * 2 + z_3 * 3
* Посчитаем число корректных предсказаний
gen z_correct = . if e(sample)
replace z_correct = 1 if ((ciga == z_class) & e(sample))
replace z_correct = 0 if ((ciga != z_class) & e(sample))
	* Доля верных предсказаний (количество единиц делить на обхем выборки)
tab z_correct if e(sample)
	* Навивный прогноз (количество индивидов, принадлежащих к самой распространенной)
	* категории следует поделить на общее число индивидов
tab ciga if e(sample)
* Проверим возможность исключить из модели переменную на образование
mlogit ciga c.age c.age#c.age i.marriage c.children c.BMI c.BMI#c.BMI work if(male == 1 & age >= 18 & educ != .)
estimates store mlogit_model_R
	* Смотрим на p-value и отвергаем нулевую гипотезу на любом разумном уровне значимости
lrtest mlogit_model_R mlogit_model_F
* Посмотрим на изменения relative risk ratio (отношение вероятнотей попадания в различные категории)
* которые показывают, как изменяется вероятность P(z=q)/P(z=базовая категория)
mlogit ciga c.age c.age#c.age i.educ i.marriage c.children c.BMI c.BMI#c.BMI work if(male == 1 & age >= 18)
mlogit, rrr
* ЗАДАНИЯ
* 1. Повторите оценивание модели включив переменную на место жительство status, а также проверив, при помощи
*    lr теста, возможность её исключения из модели.
* 2. Посчитайте отношения риска rrr по переменной age
* 3. Постройте mlogit модель для отрасли, в которой работает индивид (переменная) u_occup08,
*    предварительно выделив не более пяти категорий, в каждой из которых должно быть не менее 100 индивидов
********
* РАЗДЕЛ №3. Тестирование независимости от иррелевантных альтернатив IIA
********
* Подготовимся к проведению теста
	* Переменная на образование
generate educ_3 = 1 if educ == 3
replace educ_3 = 0 if (educ != 3) & (educ != .)
generate educ_2 = 1 if educ == 2
replace educ_2 = 0 if (educ != 2) & (educ != .)
generate educ_1 = 1 if educ == 1
replace educ_1 = 0 if (educ != 1) & (educ != .)
	* Переменная на квадрат возраста и BMI
gen age2 = age ^ 2
gen BMI2 = BMI ^ 2
	* Взаимодействие между детьми и браком
gen children_marriage = c.children * c.marriage
	* Оценим модель
mlogit ciga c.age c.age2 educ_1 educ_2 educ_3 marriage c.children c.BMI c.BMI2 work if(male == 1 & age >= 18)
* Нулевая гипотеза отвергается не для всех категорий
mlogtest, iia
