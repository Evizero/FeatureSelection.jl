ty = [21.0, 21.0, 22.8, 21.4, 18.7, 18.1, 14.3,
      24.4, 22.8, 19.2, 17.8, 16.4, 17.3, 15.2,
      10.4, 10.4, 14.7, 32.4, 30.4, 33.9, 21.5,
      15.5, 15.2, 13.3, 19.2, 27.3, 26.0, 30.4,
      15.8, 19.7, 15.0, 21.4]

tyh = [22.42049, 22.42049, 22.03658, 16.12434,
       16.66181, 13.66730, 17.12251, 20.80806,
       22.57405, 22.57405, 22.57405, 16.04756,
       16.04756, 16.04756, 14.97260, 15.51008,
       17.27607, 23.80257, 30.32907, 24.87752,
       20.88484, 13.66730, 16.66181, 21.11519,
       16.12434, 23.80257, 26.48995, 21.42232,
       24.87752, 20.27058, 19.65632, 24.03292]

tn=2

@test round(aic(ty,tyh,tn),1) == 190.8