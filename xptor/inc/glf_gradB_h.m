      real*8 gradBhp1,gradBhp3
      real*8 gradBhr11,gradBhr13,gradBhr33
      real*8 gradBhu1,gradBhu3,gradBhu33
      real*8 ave_gradBhp1(ns,nb,nb),ave_gradBhp3(ns,nb,nb)
      real*8 ave_gradBhr11(ns,nb,nb),ave_gradBhr13(ns,nb,nb)
      real*8 ave_gradBhr33(ns,nb,nb)
      real*8 ave_gradBhu1(ns,nb,nb),ave_gradBhu3(ns,nb,nb)
      real*8 ave_gradBhu33(ns,nb,nb)
      common /gradBh/
     > gradBhp1,gradBhp3,gradBhr11,gradBhr13,gradBhr33,
     > gradBhu1,gradBhu3,gradBhu33,
     > ave_gradBhp1,ave_gradBhp3,
     > ave_gradBhr11,ave_gradBhr13,ave_gradBhr33,
     > ave_gradBhu1,ave_gradBhu3,ave_gradBhu33 
  
