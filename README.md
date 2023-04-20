# primate_metabolomics

### Сложность: ⭐⭐⭐

## **Дано**

Концентрации метаболитов, измеренных в мозге, почках и мышцах человека, шимпанзе и макаки. Про каждый образец известен пол, возраст и отличается ли он по технологии — колонка *outlier* (файлы [pbio.1001871.s016](https://docs.google.com/spreadsheets/d/1FljSrt6Fc0xPZTpEvE3tadVf3rEeHFG0/edit?usp=sharing&ouid=102152803830635899757&rtpof=true&sd=true) и [pbio.1001871.s022](https://docs.google.com/spreadsheets/d/1D3rI5Pf1b2sgQl_eruNYJ3FVq4mInogk/edit?usp=share_link&ouid=102152803830635899757&rtpof=true&sd=true)).

Локальные копии таблиц сохранены в папке data проекта.

# Условные обозначения
|Название	|Определение|
| ------------- | ------------- |
|Human	|Образцы человека|
|Chimp	|Образцы шимпанзе|
|Rhesus	|Образцы макаки|
|PFC	|Префронтальная кора мозга|
|V1	|Первичная визуальная кора|
|CBC	|Мозжечок|
|Kidney	|Почки|
|Muscle	|Мышцы|
|lc_pos-Peak_	|Типичное обозначение одного метаболита|
|HMDB	|Идентификатор метаболита в БД HMDB, для тех метаболитов,<br />которые удалось четко идентифицировать|

## Что можно посмотреть

- Можно ли по метаболитам различить каждого из приматов? Во всех тканях?
- Человек-специфичные метаболиты. Есть ли они, в каких тканях и сколько? 
Сильно ли различие?
- Аналогично можно посчитать шимпанзе- и макака-специфичные метаболиты
- Есть ли влияние пола на поведение метаболитов внутри каждого вида в каждой ткани?
- Есть ли влияние возраста на поведение метаболитов внутри каждого вида в каждой ткани?

### Препроцессинг

1. Посмотрите на обе таблицы.
2. Определите, какой формат данных использован в каждом случае (широкий или длинный).
3. Попробуйте привести таблицы к одинаковому формату
4. Попробуйте соединить обе таблицы в одну.

### Результаты
![shiny app demonstration](/files/shiny_screen.png)