* Потанин Богдан Станиславович
* Микроэконометрика
* Семинар №3
* Выбор оптимальной спецификации модели бинарного выбора
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
* РАЗДЕЛ №1. Тестирование гипотезы о нормальном распределении в пробит модели
********
* Установим пакет для проведения LM теста на нормальность случайной ошибки
. ssc install skprobit
* Данный пакет не поддерживает указание типов переменных, поэтому нам придется создать новые 
	* Переменная на образование
generate educ_3 = 1 if educ == 3
replace educ_3 = 0 if (educ != 3) & (educ != .)
generate educ_2 = 1 if educ == 2
replace educ_2 = 0 if (educ != 2) & (educ != .)
generate educ_1 = 1 if educ == 1
replace educ_1 = 0 if (educ != 1) & (educ != .)
	* Переменная на квадрат возраста
gen age2 = age ^ 2
	* Взаимодействие между детьми и браком
gen children_marriage = c.children * c.marriage
* Сперва оценим пробит модель
probit work age age2 educ_1 educ_2 educ_3 children marriage children_marriage BMI  if (male == 1 & age >= 18)
* Осуществим тест и обратим внимание, что в силу малого значения p-value нулевая гипотеза
* о том, что распределение случайной ошибки является нормальным, отвергается на любом
* разумном уровне значимости.
skprobit work age age2 educ_1 educ_2 educ_3 children marriage children_marriage BMI  if (male == 1 & age >= 18)
	* ЗАДАНИЯ
	* 1. Повторите тест, добавив переменную на место жительства (status)
	* 2. Перечислите последствия, к которым может привести нарушение гипотезы о нормальном
	*    распределение случайной ошибки
* Непараметрические модели бинарного выбора, реализованные в stata, описаны в статье:
* https://journals.sagepub.com/doi/pdf/10.1177/1536867X0800800203
* Установим необходимый пакет через поиск, найдя указанную выше статью
* . ssc install snp
* Воспользуемся непараметрическим ядерным методом Klein and Spady (очень долго считает, на семинаре не запускать)
* sml work age age2 educ_1 educ_2 educ_3 children marriage children_marriage BMI  if (male == 1 & age >= 18)
* Используем полу-непараметрический подход Galland and Nychka и обратим внимание, что здесь также
* на любом разумном уровне значимости отвергается нулевая гипотеза о нормальном распределении.
* Данная функция не очень надежна и не входит в перечень оффициальных пакетов stata
* snp work age age2 educ_1 educ_2 educ_3 children marriage children_marriage BMI  if (male == 1 & age >= 18), order (3)
********
* РАЗДЕЛ №2. Тестирование гипотезы о наличии гетероскедастичности и её непосредственный учет в пробит модели
********
* Оценим модель с учетом гетероскедастичности предполагая, что дисперсия случайной ошибки разнится в зависимости
* от количества детей, возраста и брака, то есть sigma=exp(g1*age+g2*children+g3*marriage). Нулевая гипотеза об отсутствии гетероскедастичности отвергается.
hetprobit work c.age c.age#c.age i.educ c.children i.marriage i.marriage#c.children c.BMI  if (male == 1 & age >= 18), het(c.age c.children i.marriage)
* Сохраним результаты вычислений и выборку, на которой строилась модель
 estimates store probit_model
 gen probit_model_sample = e(sample)
 * Посчитаем предельные эффекты
	* Средний предельный эффект
margins, dydx(*)
	* Предельные эффекты для среднего индивида
margins, dydx(*) atmeans
	* Предельные эффекты для конкретного индивида
		* Тридцателетний индивид с двумя детьми и высшим образованием
margins, dydx(*) atmeans at(age = 30 educ = 3 children = 2 marriage = 1 BMI = 25)
		* Пятидесятилетний индивид без детей со средним образованием и несколькими средними характеристиками
margins, dydx(*) atmeans at(age = 50 educ = 1 children = 0)
* Оценим предиктивное качество модели
	* Посмотрим на число верных предсказаний
estat class, cutoff(0.5)
* Рассмотрим предсказанные вероятности
predict probit_prob, pr
sum probit_prob
histogram probit_prob
	* Предскажем вероятность того, что индивид с определенными характеристиками 
	* будет работать. Воспользуемся полученными оценками коэффициентов.
	* Для удаления ранее созданных в mata переменных используйте "mata: mata clear" без кавычек
. mata normal(                             /* считаем латентную переменную    */
			  (-7.494479 +                 /* константа                       */
			  (30 * 0.4121106) +           /* 30 лет                          */
              (30 ^ 2 *  -0.0054454) +     /* 30 лет в квадрате               */
			  0.5868622 +                  /* среднее специальное образование */
			  3 * 0.3703634 +              /* трое детей                      */
			  1.274366 +                   /* женатый                         */
			  (3 * 1) * -0.4384118 +       /* взаимодействие                  */
			  0.0256189 * 23) /            /* индекс массы тела               */
										   /* считаем дисперсию               */
			  exp((30 * .0139845 + 3 * 0.0966299  -.0987028))
			  )             
* Проверим гипотезу о наличии гетероскедастичности при помощи LM теста
	* Найдем реализации оценок латентной переменной и функции от неё
predict probit_latent, xb
gen F_probit_latent = normal(probit_latent)
gen f_probit_latent = normalden(probit_latent)
	* Посчитаем обобщенные остатки
gen generalized_residuals = ((work - F_probit_latent) / (F_probit_latent * (1 - F_probit_latent))) * f_probit_latent
	* Вычисли производные по параметрам
		* По параметрам латентной переменной (stata не поддерживает слишком длинные названия переменных)
gen derivative_beta_age = generalized_residuals * c.age
gen derivative_beta_age2 = generalized_residuals * c.age * c.age
gen derivative_beta_educ_1 = generalized_residuals * educ_1
gen derivative_beta_educ_2 = generalized_residuals * educ_2
gen derivative_beta_educ_3 = generalized_residuals * educ_3
gen derivative_beta_children = generalized_residuals * children
gen derivative_beta_marriage = generalized_residuals * marriage
gen derivative_beta_m_c = generalized_residuals * marriage * children
gen derivative_beta_BMI = generalized_residuals * BMI
	* По параметрам переменных, влияющих на дисперсию
gen derivative_alpha_age = generalized_residuals * probit_latent * age	
gen derivative_alpha_children = generalized_residuals * probit_latent * children	
gen derivative_alpha_marriage = generalized_residuals * probit_latent * marriage	
	* Создаем вектор из единиц
gen ones_vector = age / age
	* Строим регрессию на вектор из единиц без константы, в качеств независимых переменных используя частные производные
	* логарифма функция правдоподобия по оцениваемым параметрам, учитывая, что коэффициенты при переменных, влияющих
	* на дисперсию случайной ошибки, в соответствии с нулевой гипотезой равняются нулю
reg ones_vector derivative_beta_age derivative_beta_age2 derivative_beta_educ_1 derivative_beta_educ_2 derivative_beta_educ_3 derivative_beta_children derivative_beta_marriage derivative_beta_m_c derivative_beta_BMI derivative_alpha_age derivative_alpha_children derivative_alpha_marriage if (male == 1 & age >= 18), noconstant 
	* Считаем статику теста
gen LM_stat = e(r2) * e(N)
	* Вычисляем p-value теста учитывая, что число степеней свободы распределения тестовой статистики равняется количеству
	* влияющих на дисперсию случайной ошибки переменных, то есть в данном случае трем. Поскольку p-value крайне мал, то
	* нулевая гипотеза об отсутствии гетероскедастичности отвергается на любом разумном уровне значимости.
disp 1 - chi2(3, LM_stat)
	* ЗАДАНИЯ
	* 1. Сравните все полученные результаты, в частности, предельные эффекты и их стандартные ошибки, с теми, что
	*    были получены в пробит модели без учета гетероскедастичности
	* 2. Повторите оценивание модели учтя влияние лишь брака на дисперсию случайной ошибки при помощи LR и LM тестов. Провертье, значимо ли его влияние?
	* 3. Повторите оценивание модели учтя влияние BMI и квадрата возраста на дисперсию ошибки при помощи LR и LM тестов. Провертье, значимо ли их влияние?
********
* РАЗДЕЛ №3. Тестирование совместной значимости коэффициентов
********
* Проведем тест на совместную значимость переменных на образование при помощи LR теста
	* Оценим полную модель
probit work c.age c.age#c.age i.educ c.children i.marriage i.marriage#c.children c.BMI if (male == 1 & age >= 18)
estimates store m_F
gen sample_F = e(sample)
	* Оценим ограниченную модель по той же выборке, что и полную
probit work c.age c.age#c.age c.children i.marriage i.marriage#c.children c.BMI if (sample_F == 1)
estimates store m_R
	* Нулевая гипотеза LR теста отвергается, а значит предпочтительна полная модель,
	* следовательно, коэффициенты при переменных на уровень образования совместное значимы
lrtest m_F m_R
* Теперь воспользуемся Wald тестом для проверки гипотезы о значимости переменных на возраст
probit work c.age c.age#c.age i.educ c.children i.marriage i.marriage#c.children c.BMI if (male == 1 & age >= 18)
test c.age c.age#c.age
	* ЗАДАНИЯ
	* 1. Проверьте, целессобразно ли исключить из модели переменную на количество детей и брак
	* 2. Проверьте, целессобразно ли включить в модель квадрат BMI и куб возраста.
********
* РАЗДЕЛ №4. Выбор между моделями, оцененными по разными выборкам
********
* При помощи LR теста проверим, можно ли оценивать единую модель для женатых
* и холостых мужчин
	* Сперва оценим общую модель и сохраним значение логарифма правдоподобия
probit work c.age c.age#c.age i.educ c.children i.marriage i.marriage#c.children c.BMI if (male == 1 & age >= 18)
		* Логарифм правдоподобия ограниченной модели 
gen m_R = e(ll)
	* Теперь оценим полную модель, по сути оценив отдельно модели для женатых и холостых
probit work c.age c.age#c.age i.educ c.children c.BMI if (male == 1 & age >= 18 & marriage == 1)
		* Логарифм правдоподобия полной модели для женатых *
gen m_F_1 = e(ll)
probit work c.age c.age#c.age i.educ c.children c.BMI if (male == 1 & age >= 18 & marriage == 0)
		* Логарифм правдоподобия полной модели для холостых *
gen m_F_0 = e(ll)
		* Логарифм правдоподобия полной модели *
gen m_F = m_F_0 + m_F_1
	* Рассчитаем статистику теста
gen lr_stat = 2 * (m_F - m_R)
	* Найдем число степеней теста, равное числу наложенный ограниченной моделью ограничений
	* на оцениваемые параметры. Поскольку ограниченная модель предполагает, что все коэффициенты,
	* включая за исключением константы и числа детей, одинаковы для мужчин и для женщин, то число степеней свободы будет равняться
	* числу независимых переменных, за вычетом этих двух переменных, то есть 7.
	* Посчитаем p-value теста и обратим внимание, что на любом разумном уровне значимости нулевая гипотеза
	* отвергается, а значит строить общую для женатых и холостых мужчин модель нельзя.
	* ЗАДАНИЯ
	* 1. Проверьте, можно ли оценивать единую модель людей старше 35 лет и младше 35 лет
	* 2. Проверьте, можно ли оценивать единую модель для людей с детьми и без детей
disp 1 - chi2(7, lr_stat)
********
* РАЗДЕЛ №5. Оценивание параметров системы бинарных уравнений
********
* Оценим совместную модель для вероятности занятости и брака и обратим внимание, что нулевая гипотеза о равенстве 
* коэффициента корреляции между случайными ошибками (rho0) нулю отвергается, а значит модели нужно оценивать вместе.
* Также, обратите внимание на изменение знака реализации оенки коэффициента при переменной на замужество.
biprobit (work = c.age c.age#c.age i.educ c.children i.marriage i.marriage#c.children c.BMI) (marriage = i.educ c.age c.age#c.age c.BMI c.BMI#c.BMI) if (male == 1 & age >= 18)
* Сохраним информацию
estimates store biprobit_model
 gen biprobit_model_sample = e(sample)
* Предскажем различные вероятности, а именно,
* https://www.stata.com/manuals/rbiprobitpostestimation.pdf
	* Что индивид работает и в браке
predict biprobit_prob_1_1, p11
sum biprobit_prob_1_1
	* Что индивид работает и холост
predict biprobit_prob_1_0, p10
sum biprobit_prob_1_0
	* Что индивид не работает и в браке
predict biprobit_prob_0_1, p01
sum biprobit_prob_0_1
	* Что индивид не работает и холост
predict biprobit_prob_0_0, p00
sum biprobit_prob_0_0
	* Что индивид работает, при условии, что он в браке/холост
predict biprobit_prob_1_cond, pcond1
sum biprobit_prob_1_cond
sum biprobit_prob_1_cond if (male == 1 & age >= 18 & marriage == 1)
sum biprobit_prob_1_cond if (male == 1 & age >= 18 & marriage == 0)
hist(biprobit_prob_1_cond)
	* Что индивид работает
predict biprobit_prob_1, pmarg1
sum biprobit_prob_1
hist(biprobit_prob_1)
******** 
* РАЗДЕЛ №6. Выбор между моделями на основании критериев AIC и BIC 
********
* Осуществим выбор между пробит, логит, и гетпробит моделями на основании критерия AIC и BIC
	* Пробит модель
probit work c.age c.age#c.age i.educ c.children i.marriage i.marriage#c.children c.BMI if (male == 1 & age >= 18)
estimates store model_probit
	* Логит модель
logit work c.age c.age#c.age i.educ c.children i.marriage i.marriage#c.children c.BMI if (male == 1 & age >= 18)
estimates store model_logit
	* Пробит модель с учетом гетероскедастичности
hetprobit work c.age c.age#c.age i.educ c.children i.marriage i.marriage#c.children c.BMI  if (male == 1 & age >= 18), het(c.age c.children i.marriage)
estimates store model_hetprobit
	* Система бинарных уравнений
biprobit (work = c.age c.age#c.age i.educ c.children i.marriage i.marriage#c.children c.BMI) (marriage = i.educ c.age c.age#c.age c.BMI c.BMI#c.BMI) if (male == 1 & age >= 18)
estimates store model_biprobit
	* Сравним модели и учтем, что сравнивать biprobit модель с остальными в данном случае нельзя, потому что
	* она оценена, фактически, по другой выборке (подумайте, как сравнить biprobit модель с остальными).
estimates stats model_probit model_logit model_hetprobit model_biprobit
	* ЗАДАНИЯ
	* 1. Подумайте, как сравнить biprobit модель с остальными
	* 2. Сравните модели по out of sample prediction точности
