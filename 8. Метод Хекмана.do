* Потанин Богдан Станиславович
* Микроэконометрика
* Семинар №8
* Метод Хекмана
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
mata X = (rnormal(n, 1, 0, 1), rnormal(n, 1, 0, 1), rnormal(n, 1, 0, 1))
	* Коэффициенты основного уравнения
mata beta_0 = 1
mata beta_1 = 2
mata beta_2 = -2
* Коэффициенты уравнения отбора
mata gamma_0 = 1
mata gamma_1 = 1
mata gamma_2 = 1
	* Случайная ошибка
mata sigma = 2.5
mata rho = 0.8
		* Симулируем случайные ошибки из двумерного нормального распределения
		* при помощи двух стандартных нормальных случайных величин
mata z_1 = rnormal(n, 1, 0, 1)
mata z_2 = rnormal(n, 1, 0, 1)
mata epsilon_y = sigma :* z_1
mata epsilon_z = rho :* z_1 + sqrt(1 - rho ^ 2) :* z_2
	* Зависимая переменная основного уравнения
mata y_star = beta_0 :+ beta_1 :* X[,1] :+ beta_2 :* X[,2] + epsilon_y
	* Зависимая переменная уравнения отбора
mata z_star = gamma_0 :+ gamma_1 :* X[,1] :+ gamma_2 :* X[,3] + epsilon_z
	* Наблюдаемые значения зависимых переменных
	* Наблюдаемые значения зависимой переменной
mata y = y_star
mata z = z_star
. mata
for(i = 1; i <= n; i++)
{
	if(z_star[i,1] >= 0)
	{
		z[i,1] = 1
	}
	if(z_star[i,1] < 0)
	{
		z[i,1] = 0
		y[i,1] = .
	}
}
end
* Запрограммируем функцию правдоподобия
. mata
void lnL(todo, x, X, Y, Z, y, g, H)
{
	sigma_m = x[1]
	rho_m = x[2]
	
	beta_0_m = x[3]
	beta_1_m = x[4]
	beta_2_m = x[5]
	
	gamma_0_m = x[6]
	gamma_1_m = x[7]
	gamma_2_m = x[8]
	
	n = rows(X)
	
	xb_y = beta_0_m :+ beta_1_m :* X[,1] :+ beta_2_m :* X[,2]
	xb_z = gamma_0_m :+ gamma_1_m :* X[,1] :+ gamma_2_m :* X[,3]
	
	y_res = Y - xb_y
	
	L_vector = J(n, 1, .)

	for(i = 1; i <= n; i++)
	{
		if (Z[i,1] == 1)
		{
			mean_cond = (rho_m / sigma_m) * y_res[i,1]
			sd_cond = (1 - rho_m ^ 2) ^ 0.5
			L_vector[i,1] = (1 - normal((-xb_z[i,1] - mean_cond) / sd_cond)) * normalden(y_res[i,1], 0, sigma_m)
		}
		if (Z[i,1] == 0)
		{
			L_vector[i,1] = normal(-xb_z[i,1])
		}
	}
	
	y = sum(log(L_vector))
}
end
* Запустим оптимизатор
* Возьмем начальные точки близкие к истинным, обычно они берутся из
* двухшаговой процедуры
. mata
S = optimize_init()
optimize_init_evaluator(S, &lnL())
optimize_init_params(S, (2.5, 0.8, 1.1, 2.1, -2.1, 1.1, 1.1, 1.1))
optimize_init_argument(S, 1, X)
optimize_init_argument(S, 2, y)
optimize_init_argument(S, 3, z)
optimize_init_conv_maxiter(S, 1000)
x = optimize(S)
H = optimize_result_Hessian(S)
end
* Посмотрим на полученные значения и сравним их с истинными
mata (x \ (sigma, rho, beta_0, beta_1, beta_2, gamma_0, gamma_1, gamma_2))
* ЗАДАНИЯ
* 1. Сравните полученные оценки с полученными при помощи МНК
* 2. Посмотрите, как изменится точность оценок если не будут соблюдены exclusion restrictions.
*    Например, сделайте так, чтобы и основное уравнение и уравнение отбора включали лишь X1 и X2
*    в качестве независимых переменных.
* 3. Посмотрите, как изменится точность оценок при нарушении допущения о совместном нормальном
*    распределении случайных ошибок
********
* РАЗДЕЛ №2. Приложение метода Хекмана к реальным данным
********
* https://www.stata.com/manuals13/rtobitpostestimation.pdf
* Добавим переменную на зарплату
gen wage = uj13_2
replace wage = . if uj13_2 > 100000000
replace wage = . if work == 0
gen wage_ln = log(wage + 1)
* И на стаж
gen seniority = uj161_3y
replace seniority = . if seniority > 1000
* Оценим обыкновенную регрессию
reg wage_ln c.seniority c.seniority#c.seniority i.educ if((male == 1) & (age >= 18))
* Осуществим двухшаговую процедуру Хекмана
heckman wage_ln c.seniority c.seniority#c.seniority i.educ if((male == 1) & (age >= 18)), select(work = c.age c.age#c.age i.educ i.marriage c.children c.children#c.children) twostep
* Воспользуемся методом Хекмана, основанном на методе максимального правдоподобия
heckman wage_ln c.seniority c.seniority#c.seniority i.educ if((male == 1) & (age >= 18)), select(work = c.age c.age#c.age i.educ i.marriage c.children c.children#c.children)
 * Посчитаем предельные эффекты
	* На безусловное математическое ожидание E(y_star)
margins, dydx(*)
	* На условное математическое ожидание E(y_star|z=1)
margins, dydx(*) predict(e(0,.))
* Предскажем
	* Значения y_star
predict wage_xb, xb
	* Значение y_star|z=1
predict wage_cond, ycond
	* Значение y с присвоением 0 для незанятых
predict wage_cond, yexpected
* ЗАДАНИЯ
* 1. Проверьте возможность исключить из модели переменные на образование
* 2. Сравните предсказательную out of sample точность модели с образованием и без него
********
* РАЗДЕЛ №3. Модель Хекмана для порядковой пробит модели
********
* Создадим переменную на удовлетворенность жизнью
* Уровни удовлетворенность работой
	* Полностью удовлетворены
generate worksat = 5 if uj1_1_1 == 1
	* Скорее удовлетворены
replace worksat = 4 if uj1_1_1 == 2
	* И да, и нет
replace worksat = 3 if uj1_1_1 == 3
	* Не очень удовлетворены
replace worksat = 2 if uj1_1_1 == 4
	* Совсем не удовлетворены
replace worksat = 1 if uj1_1_1 == 5
* Воспользуемся методом Хекмана
heckoprobit worksat c.seniority c.seniority#c.seniority i.educ if((male == 1) & (age >= 18)), select(work = c.age c.age#c.age i.educ i.marriage c.children c.children#c.children)
 * Посчитаем предельные эффекты
	* На безусловную вероятность того или иного уровня удовлетворенности жизнью
margins, dydx(*)
* ЗАДАНИЯ
* 1. Проверьте, можно ли оценивать единую модель для мужчин и для женщин
