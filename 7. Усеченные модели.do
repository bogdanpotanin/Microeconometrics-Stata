* Потанин Богдан Станиславович
* Микроэконометрика
* Семинар №7
* Усеченные модели
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
* РАЗДЕЛ №1. Симуляционный анализ модели тобита
********
mata: mata clear
* Приступим к симуляциям
	* Количество наблюдений
mata n = 10000
	* Матрица независимых переменных
mata X = (runiform(n, 1, 1, 2), runiform(n, 1, 1, 2))
	* Коэффициенты
mata beta_0 = 1
mata beta_1 = 2
mata beta_2 = -2
	* Случайная ошибка
mata sigma = 2.5
mata epsilon = rnormal(n, 1, 0, sigma)
	* Зависимая переменная
mata y_star = beta_0 :+ beta_1 :* X[,1] + beta_2 :* X[,2] + epsilon
	* Усечения
mata tr_low = -1
mata tr_up = 1.5
	* Наблюдаемые значения зависимой переменной
mata y = y_star
. mata
for(i = 1; i <= n; i++)
{
	if((y_star[i,1] < tr_low))
	{
		y[i,1] = tr_low
	}
	if((y_star[i,1] > tr_up))
	{
		y[i,1] = tr_up
	}
}
end
* Запрограммируем функцию правдоподобия
. mata
void lnL(todo, x, X, Y, tr_low, tr_up, y, g, H)
{
	beta_0_m = x[1]
	beta_1_m = x[2]
	beta_2_m = x[3]
	sigma_m = x[4]
	
	n = rows(X)
	
	xb = beta_0_m :+ beta_1_m :* X[,1] + beta_2_m :* X[,2]
	
	L_vector = J(n, 1, .)
	
	for(i = 1; i <= n; i++)
	{
		if ((Y[i,1] == tr_low))
		{
			L_vector[i,1] = normal((tr_low :- xb[i,1]) :/ sigma_m)
		}
		if ((Y[i,1] == tr_up))
		{
			L_vector[i,1] = 1 - normal((tr_up :- xb[i,1]) :/ sigma_m)
		}
		if (((Y[i,1] > tr_low) & (Y[i,1] < tr_up)))
		{
			L_vector[i,1] = normalden(Y[i,1] :- xb[i,1], 0, sigma_m)
		}
	}
	
	y = sum(log(L_vector))
}
end
* Запустим оптимизатор
. mata
S = optimize_init()
optimize_init_evaluator(S, &lnL())
optimize_init_params(S, (0, 0, 0, 1))
optimize_init_argument(S, 1, X)
optimize_init_argument(S, 2, y)
optimize_init_argument(S, 3, tr_low)
optimize_init_argument(S, 4, tr_up)
optimize_init_conv_maxiter(S, 1000)
x = optimize(S)
H = optimize_result_Hessian(S)
end
* Посмотрим на полученные значения и сравним их с истинными
mata (x \ (beta_0, beta_1, beta_2, sigma))
* Теперь посмотрим на усеченную регрессию
. mata
void lnL_new(todo, x, X, Y, tr_low, tr_up, y, g, H)
{
	beta_0_m = x[1]
	beta_1_m = x[2]
	beta_2_m = x[3]
	sigma_m = x[4]
	
	n = rows(X)
	
	xb = beta_0_m :+ beta_1_m :* X[,1] + beta_2_m :* X[,2]
	
	L_vector = J(n, 1, 1)
	
	for(i = 1; i <= n; i++)
	{
		if (((Y[i,1] > tr_low) & (Y[i,1] < tr_up)))
		{
			L_vector[i,1] = normalden(Y[i,1] :- xb[i,1], 0, sigma_m) :/ (normal((tr_up :- xb[i,1]) :/ sigma_m) :- normal((tr_low :- xb[i,1]) :/ sigma_m))
		}
	}
	
	y = sum(log(L_vector))
}
end
* Запустим оптимизатор
. mata
S = optimize_init()
optimize_init_evaluator(S, &lnL_new())
optimize_init_params(S, (0, 0, 0, 1))
optimize_init_argument(S, 1, X)
optimize_init_argument(S, 2, y)
optimize_init_argument(S, 3, tr_low)
optimize_init_argument(S, 4, tr_up)
optimize_init_conv_maxiter(S, 1000)
x = optimize(S)
H = optimize_result_Hessian(S)
end
* Посмотрим на полученные значения и сравним их с истинными
mata (x \ (beta_0, beta_1, beta_2, sigma))
********
* РАЗДЕЛ №2. Приложение модели Тобина к реальным данным
********
* https://www.stata.com/manuals13/rtobitpostestimation.pdf
* Добавим переменную на зарплату
gen wage = uj13_2
replace wage = . if uj13_2 > 100000000
replace wage = 0 if work == 0
gen wage_ln = log(wage + 1)
* Оценим обыкновенную регрессию
reg wage_ln c.age c.age#c.age i.educ c.BMI if((male == 1) & (age >= 18))
* Посчитаем модель Тобина для совершеннолетних мужчин
tobit wage_ln c.age c.age#c.age i.educ c.BMI if((male == 1) & (age >= 18)), ll(0)
 * Посчитаем предельные эффекты
	* На логарфим зарплаты (y_star)
margins, dydx(*)
	* На условное математическое ожидание (y=y_star|tr_low<=y<=tr_up)
margins, dydx(*) predict(e(0,.))
	* На вероятность не попасть под усечение (p(y_star>0))
margins, dydx(*) predict(p(0,.))
* Предскажем
	* Значения y_star
predict wage_xb, xb
	* Вероятность того, что не произойдет усечение
predict wage_pr, pr(0,9999999999)
	* Значение y (хоть название и ystar но предсказывает игрек без звездочки)
predict wage_y, ystar(0,9999999999)
* Оценим предиттивную силу модели посмотрев 
	* на корреляцию между истинными и предсказанными значениями среди индивидов с ненулевой заработной платой
correlate wage_xb wage_ln if (e(sample) & (wage_ln > 0 ))
	* Сравним с МНК
reg wage_ln c.age c.age#c.age i.educ c.BMI if((male == 1) & (age >= 18))
predict wage_reg, xb
correlate wage_reg wage_ln if (e(sample) & (wage_ln > 0 ))
	* Однако, конечно, для адекватного сравнения следует использовать out of sample prediction
* Сравним точность моделей используя предсказания за пределами выборки
	* Сгенерируем случайную величину из равногомерного распределения
gen train_val = runiform()
	* Отберем около 66% выборки в качестве тренировочной
gen train_ind = train_val > 0.66
	* Оценим Тобит модель по тренировочной выборке и посмотрим на mse
tobit wage_ln c.age c.age#c.age i.educ c.BMI if((male == 1) & (age >= 18) & (train_ind == 1)), ll(0)
predict pred_tobit, ystar(0,9999999999)
gen se_tobit = (pred_tobit - wage_ln) ^ 2
mean se_tobit if ((wage_ln > 0)  & (!train_ind) & (male == 1) & (age >= 18))
mean se_tobit if ((!train_ind) & (male == 1) & (age >= 18))
twoway scatter pred_tobit wage_ln if ((wage_ln > 0)  & (!train_ind) & (male == 1) & (age >= 18))
	* Оценим обычную модель по тренировочной выборке и посмотрим на корреляцию
reg wage_ln c.age c.age#c.age i.educ c.BMI if((male == 1) & (age >= 18) & (train_ind == 1))
predict pred_reg, xb
gen se_reg = (pred_reg - wage_ln) ^ 2
mean se_reg if ((wage_ln > 0)  & (!train_ind) & (male == 1) & (age >= 18))
mean se_reg if ((!train_ind) & (male == 1) & (age >= 18))
* ЗАДАНИЯ
* 1. Проверьте возможность исключить из модели переменные на образование
* 2. Сравните предсказательную out of sample точность модели с образованием и без него
********
* РАЗДЕЛ №3. Усеченная регрессия
********
truncreg wage_ln c.age c.age#c.age i.educ c.BMI if((male == 1) & (age >= 18)), ll(0)* Посчитаем предельные эффекты
	* На логарфим зарплаты (y_star)
margins, dydx(*)
	* На условное математическое ожидание (y=y_star|tr_low<=y<=tr_up)
margins, dydx(*) predict(e(0,.))
	* На вероятность не попасть под усечение (p(y_star>0))
margins, dydx(*) predict(p(0,.))
* ЗАДАНИЯ
* 1. Вычислите mse по out of sample prediction
* 2. При помощи mata напишите собственную функцию для усеченной регрессии
