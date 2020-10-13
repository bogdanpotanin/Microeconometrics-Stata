* Потанин Богдан Станиславович
* Микроэконометрика
* Семинар №4
* Модели порядкового выбора
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
mata X = (rnormal(n, 1, 1, 2), rnormal(n, 1, -1, 3))
	* Случайная ошибка
mata epsilon = rnormal(n, 1, 0, 1)
	* Коэффициенты
mata beta_1 = 2
mata beta_2 = 3
	* Создаем латентную переменную
mata z_star = (beta_1 :* X[,1]) + (beta_2 :* X[,2]) :+ epsilon
	* Выберем порогвые значения
mata threshold_1 = -3
mata threshold_2 = 5
	* Создадим наблюдаемые значения зависимой переменной
mata z = J(n,1, .)
. mata 
for(i = 1; i <= n; i++)
{
	if(z_star[i, 1] <= threshold_1)
	{
		z[i,1] = 1
	}
	if((threshold_1 < z_star[i, 1]) & (z_star[i, 1] <= threshold_2))
	{
		z[i,1] = 2
	}
	if(z_star[i, 1] > threshold_2)
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
	threshold_1_m = x[1]
	threshold_2_m = x[2]
	
	beta_1_m = x[3]
	beta_2_m = x[4]
	
	n_m = rows(X)
	
	z_star_determ = (beta_1_m :* X[,1]) + (beta_2_m :* X[,2])
	
	lnL_vector = J(n_m, 1, .)
	
	for(i = 1; i <= n_m; i++)
	{
		if(z[i,1] == 1)
		{
			lnL_vector[i, 1] = normal(threshold_1_m - z_star_determ[i, 1])
		}
		if(z[i,1] == 2)
		{
			lnL_vector[i, 1] = normal(threshold_2_m - z_star_determ[i, 1]) - normal(threshold_1_m - z_star_determ[i, 1])
		}
		if(z[i,1] == 3)
		{
			lnL_vector[i, 1] =  1 - normal(threshold_2_m - z_star_determ[i, 1])
		}
	}
	
	y = sum(log(lnL_vector))
}
end
* Запустим оптимизатор
. mata
S = optimize_init()
optimize_init_evaluator(S, &lnL())
optimize_init_params(S, (-2, 3, 1.3, 1.8))
optimize_init_argument(S, 1, X)
optimize_init_argument(S, 2, z)
optimize_init_conv_maxiter(S, 1000)
x = optimize(S)
H = optimize_result_Hessian(S)
end
* Посмотрим на полученные значения и сравним их с истинными
mata (x \ (threshold_1, threshold_2, beta_1, beta_2))
* Убедимся, что МНК оценки окажутся смещенными
mata luinv(X'*X)*X'*z
* Но они были бы состоятельными, будь нам доступна латентная переменная
mata luinv(X'*X)*X'*z_star
* ЗАДАНИЯ
* 1. Проверьте, что случится, если вы неверное специфицируете количество категорий, например,
*    объединив первую и вторую в одну
* 2. Замените распределение случайной ошибки с нормального на логистическое и посмотрите, как
*    изменятся результаты оценивания параметров модели.
********
* РАЗДЕЛ №2. Порядкова пробит модель
********
* Создадим переменную на удовлетворенность жизнью
* Уровни удовлетворенности жизнью
	* Полностью удовлетворены
generate lifesat = 5 if uj65 == 1
	* Скорее удовлетворены
replace lifesat = 4 if uj65 == 2
	* И да, и нет
replace lifesat = 3 if uj65 == 3
	* Не очень удовлетворены
replace lifesat = 2 if uj65 == 4
	* Совсем не удовлетворены
replace lifesat = 1 if uj65 == 5
* Добавим переменную на зарплату
gen wage = uj13_2
replace wage = . if uj13_2 > 10000
replace wage = 0 if work == 0
* Убедимся, что у нас нет слишком маленьких категорий 
tab lifesat if ((male == 1) & (age >= 18))
* Построим порядковую модель для совершеннолетних мужчин
oprobit lifesat c.age c.age#c.age i.educ c.BMI c.BMI#c.BMI c.children i.marriage i.work c.wage if ((male == 1) & (age >= 18))
	* Сохраним результаты вычислений и выборку, на которой строилась модель
estimates store oprobit_model
gen oprobit_model_sample = e(sample)
 * Посчитаем предельные эффекты для каждой категории
	* Средний предельный эффект
margins, dydx(*)
	* Предельные эффекты для среднего индивида
margins, dydx(*) atmeans
	* Предельные эффекты для конкретного индивида
margins, dydx(*) atmeans at(age = 50 educ = 3 BMI = 23 children = 2 marriage = 1 work = 1 wage = 50000)
* Оценим точность модели
* Получим предсказанные вероятности
predict p0 p1 p2 p3 p4, pr
* Получим бинарные переменные на принадлежность к каждому из классов
gen z_0 = 0 if e(sample)
gen z_1 = 0 if e(sample)
gen z_2 = 0 if e(sample)
gen z_3 = 0 if e(sample)
gen z_4 = 0 if e(sample)
replace z_0 = 1 if ((p0 >= p1) & (p0 >= p2) & (p0 >= p3) & (p0 >= p4) & e(sample))
replace z_1 = 1 if ((p1 >= p1) & (p1 >= p2) & (p1 >= p3) & (p1 >= p4) & e(sample))
replace z_2 = 1 if ((p2 >= p1) & (p2 >= p2) & (p2 >= p3) & (p2 >= p4) & e(sample))
replace z_3 = 1 if ((p3 >= p1) & (p3 >= p2) & (p3 >= p3) & (p3 >= p4) & e(sample))
replace z_4 = 1 if ((p4 >= p1) & (p4 >= p2) & (p4 >= p3) & (p4 >= p4) & e(sample))
* Создадим переменную на предсказанный класс
gen z_class = z_0 * 1 + z_1 * 2 + z_2 * 3 + z_3 * 4 + z_4 * 5
* Посчитаем число корректных предсказаний
gen z_correct = . if (oprobit_model_sample)
replace z_correct = 1 if ((lifesat == z_class) & (oprobit_model_sample))
replace z_correct = 0 if ((lifesat != z_class) & (oprobit_model_sample))
	* Доля верных предсказаний (количество единиц делить на обхем выборки)
tab z_correct if e(sample)
	* Навивный прогноз (количество индивидов, принадлежащих к самой распространенной)
	* категории следует поделить на общее число индивидов
tab lifesat if e(sample)
* Проверим возможность исключить из модели переменную на образование
oprobit lifesat c.age c.age#c.age c.BMI c.BMI#c.BMI c.children i.marriage i.work c.wage if ((male == 1) & (age >= 18) & oprobit_model_sample)
estimates store oprobit_model_R
	* Смотрим на p-value и отвергаем нулевую гипотезу на уровне значимости 5%
lrtest oprobit_model_R oprobit_model
* ЗАДАНИЯ
* 1. Повторите оценивание модели включив переменную на место жительство status, а также проверив, при помощи
*    lr теста, возможность её исключения из модели.
* 2. Объедините некоторые категории так, чтобы у вас остались лишь три. Оцените полученную модель.
*    Подумайте, можно ли как-то её сравнить с предыдущей.
********
* РАЗДЕЛ №2. Порядкова логит модель и обобщенная порядковая логит модель
********
* Построим порядковую модель для совершеннолетних мужчин
ologit lifesat c.age c.age#c.age i.educ c.BMI c.BMI#c.BMI c.children i.marriage i.work c.wage if ((male == 1) & (age >= 18))
estimates store ologit_model
* Добавим отношения шансов
* При интерпретации отношений шансов учитывайте, что для любого номера категории k изменение в отношениях шансов p(z>k)/p(z<=k) остается прежним
ologit lifesat c.age c.age#c.age i.educ c.BMI c.BMI#c.BMI c.children i.marriage i.work c.wage if ((male == 1) & (age >= 18)), or
* Протестируем parallel regression lines assumption
* Установим необходимый пакет
	* Для теста бранта
search spost13_ado
	* Для обобщенной ологит модели (ищем st0097_1)
search gologit2
* Оценим обобщенную модель без parallel lines assumption ограничения
gologit2 lifesat c.age c.age#c.age i.educ c.BMI c.BMI#c.BMI c.children i.marriage i.work c.wage if ((male == 1) & (age >= 18)),  lrforce store(unconstrained)
* Вновь оценим с ограничением
gologit2 lifesat c.age c.age#c.age i.educ c.BMI c.BMI#c.BMI c.children i.marriage i.work c.wage if ((male == 1) & (age >= 18)),  lrforce pl store(constrained)
* Осуществим LR тест и убедимся, что нулевая гипотеза о соблюдении paralell lines assumption не отвергается
lrtest constrained unconstrained
* Для этой команды доступы и все прочие postestimates по аналогии с обычным пробитом и логитом
* Посмотрим на предельный эффект для индивида
margins, dydx(*) atmeans at(age = 50 educ = 3 BMI = 23 children = 2 marriage = 1 work = 1 wage = 50000)
* Проведем тест Бранта, преварительно перекодировав дамми переменные, нелинейные переменные и переменные взаимодействия
gen age2 = age ^ 2
gen BMI2 = BMI^2
generate educ_3 = 1 if educ == 3
replace educ_3 = 0 if (educ != 3) & (educ != .)
generate educ_2 = 1 if educ == 2
replace educ_2 = 0 if (educ != 2) & (educ != .)
generate educ_1 = 1 if educ == 1
replace educ_1 = 0 if (educ != 1) & (educ != .)
omodel logit lifesat age age2 educ_1 educ_2 educ_3 BMI BMI2 children marriage work wage if ((male == 1) & (age >= 18))
	* Убедимся, что нулевая гипотеза о соблюдении paralell lines assumption не отвергается
brant
