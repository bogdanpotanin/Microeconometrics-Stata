* Потанин Богдан Станиславович
* Микроэконометрика
* Семинар №6
* Моделированние частоты и вложенный выбор
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
* РАЗДЕЛ №1. Симуляционный анализ пуассоновской регрессии
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
	* В данной модели нет случайной ошибки, есть только
	* параметр пуассоновского распределения
mata lambda = exp(beta_0 :+ beta_1 :* X[, 1] + beta_2 :* X[, 2])
	* Сгенерируем зависимую переменную
mata z = J(n,1, .)
. mata 
for(i = 1; i <= n; i++)
{
	z[i, 1] = rpoisson(1, 1, lambda[i, 1])
}
end
	* Посмотрим на полученные значения
mata z
* Напишем функцию правдоподобия в соответствии с представленным выше
* процессом генерации данных
. mata
void lnL(todo, x, X, z, y, g, H)
{
	beta_0_m = x[1]
	beta_1_m = x[2]
	beta_2_m = x[3]
	
	lambda_m = exp(beta_0_m :+ beta_1_m :* X[, 1] + beta_2_m :* X[, 2])
	
	lnL_vector = exp(-lambda_m) :* (lambda_m :^ z) :/ factorial(z)
	
	y = sum(log(lnL_vector))
}
* Запустим оптимизатор
. mata
S = optimize_init()
optimize_init_evaluator(S, &lnL())
optimize_init_params(S, (0.25, -1.3, 3.5))
optimize_init_argument(S, 1, X)
optimize_init_argument(S, 2, z)
optimize_init_conv_maxiter(S, 1000)
x = optimize(S)
H = optimize_result_Hessian(S)
end
* Посмотрим на полученные значения и сравним их с истинными
mata (x \ (beta_0, beta_1, beta_2))
********
* РАЗДЕЛ №2. Применение пуассоновской регрессии для предсказания частоты езды на велосипеде
********
* Создадим переменную на езду на велосипеде
gen bycicle = um11314c
replace bycicle = . if (bycicle > 1000)
tab bycicle
* Оценим регрессию
poisson bycicle c.age c.age#c.age i.educ c.children c.marriage c.BMI i.work if(male == 1 & age >= 18)
	*Сохраним результаты вычислений и выборку, на которой строилась модель
estimates store poisson_model
gen poisson_model_sample = e(sample)
* Пуассоновская регрессия редко встречается на практике, поскольку опирается на допущение о равенстве математического ожидания
* и дисперсии зависимой переменной. Часто эту проблему именуют overdispersion, поскольку на самом деле дисперсия, как правило,
* оказывается больше среднего значения.
* Проверим, подходит ли пуассоновская регрессия в данном случае, протестировав гипотезу о том, что сумма стандартизированных
* остатков подчиняется хи-квадрат распределению с числом степеней свободы, равным разнице между числом наблюдений и регрессоров,
* включая константу. Поскольку нулевая гипотеза отвергается на любом разумно уровне значимости то, допущение о том, что пуассоновская 
* регрессия подходит - отвергается.
estat gof
 * Посчитаем предельные эффекты на lambda (среднее число событий)
	* Средний предельный эффект
margins, dydx(*)
	* Предельные эффекты для среднего индивида
margins, dydx(*) atmeans
	* Предельные эффекты для конкретного индивида
margins, dydx(*) atmeans at(age = 50 educ = 3 BMI = 23 children = 2 marriage = 1 work = 1)
* Оценим точность модели
* Получим предсказанные средние значения
predict bycicle_poisson_hat
list bycicle_poisson_hat
* Оценим точность предсказаний по mse
gen poisson_mse = (bycicle_poisson_hat - bycicle) ^ 2 if e(sample)
sum poisson_mse
* ЗАДАНИЯ
* 1. Проверьте, можем ли мы исключить переменную на высшее образование
* 2. Посчитайте предельный эффект работы на вероятность покататься хотя бы раз на велосипеде ни разу для индивида
*    с вашими характеристиками
* 3. Посчитайте предельный эффект BMI на вероятность покататься хотя бы раз на велосипеде ни разу для индивида
*    с вашими характеристиками
********
* РАЗДЕЛ №3. Симуляционный анализ отрицательной биномиальной регрессии
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
	* Случайная ошибка для lambda
mata alpha = 0.5
mata epsilon = log(rgamma(n, 1, 1 / alpha, alpha))
	* Теперь лямбда включает случайную ошибку
mata lambda = exp(beta_0 :+ beta_1 :* X[, 1] + beta_2 :* X[, 2] :+ epsilon)
	* Сгенерируем зависимую переменную
mata z = J(n,1, .)
. mata 
for(i = 1; i <= n; i++)
{
	z[i, 1] = rpoisson(1, 1, lambda[i, 1])
}
end
	* Посмотрим на полученные значения
mata z
* Напишем функцию правдоподобия в соответствии с представленным выше
* процессом генерации данных
. mata
void lnL_new(todo, x, X, z, y, g, H)
{
	beta_0_m = x[1]
	beta_1_m = x[2]
	beta_2_m = x[3]
	alpha_m = x[4]
	
	mu_m = exp(beta_0_m :+ beta_1_m :* X[, 1] + beta_2_m :* X[, 2])
	m_m = 1 :/ alpha_m
	p_m = 1 :/ (1 :+ (alpha_m :* mu_m))

	lnL_vector = log(gamma(m_m :+ z))
	lnL_vector = lnL_vector :- log(gamma(z :+ 1)) 
	lnL_vector = lnL_vector :- log(gamma(m_m)) 
	lnL_vector = lnL_vector :+ (m_m :* log(p_m)) 
	lnL_vector = lnL_vector :+ (z :* log(1 :- p_m))
	
	y = sum(lnL_vector)
}
end
* Запустим оптимизатор
. mata
S = optimize_init()
optimize_init_evaluator(S, &lnL_new())
optimize_init_params(S, (0.1, 0.1, 0.1, 0.1))
optimize_init_argument(S, 1, X)
optimize_init_argument(S, 2, z)
optimize_init_conv_maxiter(S, 100)
x = optimize(S)
H = optimize_result_Hessian(S)
end
* Посмотрим на полученные значения и сравним их с истинными
mata (x \ (beta_0, beta_1, beta_2, alpha))
********
* РАЗДЕЛ №4. Отрицательная биномиальная регрессия
********
* Обтатите внимание, что при alpha = 0 мы получаем Пуассоновскую регрессию
* Если нулевая гипотеза LR теста с H0:alpha=0 отвергается, то предпочтение отдается отрицательной биномиальной регрессии перед пуассоновской
nbreg bycicle c.age c.age#c.age i.educ c.children c.marriage c.BMI i.work if(male == 1 & age >= 18)
	*Сохраним результаты вычислений и выборку, на которой строилась модель
estimates store nbreg_model
gen nbreg_model_sample = e(sample)
 * Посчитаем предельные эффекты на lambda (среднее число событий)
 * ожидаемом количестве событий
	* Средний предельный эффект
margins, dydx(*)
	* Предельные эффекты для среднего индивида
margins, dydx(*) atmeans
	* Предельные эффекты для конкретного индивида
margins, dydx(*) atmeans at(age = 50 educ = 3 BMI = 23 children = 2 marriage = 1 work = 1)
* Оценим точность модели
* Получим предсказанные средние значения
predict bycicle_nbreg_hat
list bycicle_nbreg_hat
* Оценим точность предсказаний по mse
gen nbreg_mse = (bycicle_nbreg_hat - bycicle) ^ 2 if e(sample)
* Внутривыборочное предсказаний оказалось точней у модели Пуассона
sum nbreg_mse
sum poisson_mse
* Параметр alpha может быть и эндогенным (аналог хетпробита для отрицательной биномиальной регрессии)
gnbreg bycicle c.age c.age#c.age i.educ c.children c.marriage c.BMI i.work if(male == 1 & age >= 18), lnalpha(c.age c.BMI)
 * Посчитаем предельные эффекты на lambda (среднее число событий)
	* Средний предельный эффект
margins, dydx(*)
	* Предельные эффекты для среднего индивида
margins, dydx(*) atmeans
	* Предельные эффекты для конкретного индивида
margins, dydx(*) atmeans at(age = 50 educ = 3 BMI = 23 children = 2 marriage = 1 work = 1)
* Оценим точность модели
* Получим предсказанные средние значения
predict bycicle_gnbreg_hat
list bycicle_gnbreg_hat
* Оценим точность предсказаний по mse
gen gnbreg_mse = (bycicle_gnbreg_hat - bycicle) ^ 2 if e(sample)
sum gnbreg_mse
sum nbreg_mse
sum poisson_mse
* ЗАДАНИЯ
* 1. Проверьте для каждой из ранее оцененных моделей, можем ли мы исключить переменную на высшее образование
* 2. Сравните оцененные модели при помощи критерия AIC
* 3. Определите, какая из моделей лучше всех предсказывает вне выборки
********
* РАЗДЕЛ №5. Вложенная мультиномиальная логит модель (nested multinomial logit)
********
* На основе мануала https://www.stata.com/manuals13/rnlogit.pdf
* Обратим внимание, что (1-tau_{t})^2 коэффициент корреляции между случайными ошибками в гнезде t, правда,
* иногда возникает проблема, когда некоторые tau_{t} > 1.
* Поэтому если во всех гнездах tau_{t}=0, то мы получаем обыкновенную мультиномиал модель.
* На первом уровне у всех разные значеиня зависимых переменных, но одни и те же коэффициенты. На втором
* уровне разнятся коэффициенты, но значения переменных одинаковые.
* Создадим переменную для определения множества альтернатив
* Создадим переменную на предпочитаемый тип сигарет
use http://www.stata-press.com/data/r13/restaurant, clear
describe
qui nlogitgen type = restaurant(fast: Freebirds | MamasPizza, family: CafeEccell | LosNortenos| WingsNmore, fancy: Christophers | MadCows)
nlogittree restaurant type, choice(chosen)
nlogit chosen cost rating distance || type: income kids, base(family) || restaurant:, noconstant case(family_id)
* Автоматически предельные эффекты для данной модели stata, к сожалению, не считает
* Предскажем вероятности посещения ресторанов
	* Вероятность посещения конкретного ресторана p1 и вероятность посещения ресторана данного типа p2
predict p*, pr
	* Вероятность посещения ресторана при условии, что был выбран ресторан данного типа
predict condp, condp hlevel(2)
sort family_id type restaurant
list restaurant type chosen p2 p1 condp in 1/14, sepby(family_id) divider
