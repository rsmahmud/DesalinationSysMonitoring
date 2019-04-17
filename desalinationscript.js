// Double-Rankine Cycle OTEC Simulator
    // S. Goto 14 March, 2017 Static Property Ver.2.5 古川君の方法を使ったが，温海水流量と冷海水流量に応じた条件分岐と，温冷海水の出口温度の初期値に収束が依存するためこれで断念．

    // Configuration of OTEC plant

    var Ae1 = 87.4;   //Heat transfer area of evaporator Ae (UnitA) [m2]
    var Ae2 = 87.4;   //Heat transfer area of evaporator Ae (UnitB) [m2]
    var Ac1 = 87.4;   //Heat transfer area of condenser Ac (UnitA) [m2]
    var Ac2 = 87.4;   //Heat transfer area of condenser Ac (UnitB) [m2]

    var sigmae1 = 0.0028;  // Depth of blade of evaporator (UnitA) [m]
    var sigmae2 = 0.0028;  // Depth of blade of evaporator (UnitB) [m]
    var sigmac1 = 0.0028;  // Depth of blade of condenser (UnitA) [m]
    var sigmac2 = 0.0028;  // Depth of blade of condenser (UnitB) [m]

    var Ue1 = 2.0; // Overall heat transfer coefficient U of evaporator (UnitA) [kW/(m2 K)]
    var Ue2 = 2.0; // Overall heat transfer coefficient U of evaporator (UnitB) [kW/(m2 K)]
    var Uc1 = 2.0; // Overall heat transfer coefficient U of condenser (UnitA) [kW/(m2 K)]
    var Uc2 = 2.0; // Overall heat transfer coefficient U of condenser (UnitB) [kW/(m2 K)]

    var etaT=0.75; //タービン効率
    var etaP=0.80; //ポンプ効率
    var etaG=0.80; //発電効率

    UAe1 = (Ue1 * Ae1) * 1000;   //蒸発器伝熱性能(Unit A)[W/K]
    UAe2 = (Ue2 * Ae2) * 1000;   //蒸発器伝熱性能(Unit B)[W/K]
    UAc1 = (Uc1 * Ac1) * 1000;   //凝縮器伝熱性能(Unit A)[W/K]
    UAc2 = (Uc2 * Ae2) * 1000;   //凝縮器伝熱性能(Unit B)[W/K]


    // Physical constants
    var Cp = 4200.0; //比熱[J/(kg K)]
    // alert("Cp: "+Cp);
    var rhoa = 601.13; // Density of ammonia ρa [kg/m3]
    var ca = 4.800; // Specific heat of ammonia ca [kJ/(kg K)]
    var rhow = 1000.0; // Density of water ρw [kg/m3]
    var cw = 4.179; // Specific heat of water cw [kJ/(kg K)]

    // Simulation parameters
    var DATANUMBER = 30; // Number of Data Points for simulation
    var POINTNUMBER = 8;    // Number of Points for Double-Rankine Cycle 


    var i = 0, j = 0, k = 0, l = 0;
    //各ポイントの状態量 T:温度, P:圧力, V:比体積, h:比エンタルピー, S:比エントロピー
    var T=new Array(POINTNUMBER);
    var P=new Array(POINTNUMBER);
    var V=new Array(POINTNUMBER);
    var h=new Array(POINTNUMBER);
    var S=new Array(POINTNUMBER); 

    for(i = 0; i < POINTNUMBER; i++){
        T[i]=0;
        P[i]=0;
        V[i]=0;
        h[i]=0;
        S[i]=0;
    }

    //シミュレータへの入力
    var Twsi1 = 29.0;  //蒸発器入口温度(Unit A)
    var Twsi2 = 0.0;  //蒸発器入口温度(Unit B)
    var Twso1 = 0.0;  //蒸発器出口温度(Unit A)
    var Twso2 = 0.0;  //蒸発器出口温度(Unit B)
    var Mws = 39.72;   //温海水流量
    var Tcsi1 = 0.0;  //凝縮器入口温度(Unit A)
    var Tcsi2 = 9.0;  //凝縮器入口温度(Unit B)
    var Tcso1 = 0.0;  //凝縮器出口温度(Unit A)
    var Tcso2 = 0.0;  //凝縮器出口温度(Unit B)
    var Mcs = 38.89;   //冷海水流量
    var Mwf1 = 0.35; //作動流体流量(Unit A)[kg/s]
    var Mwf2 = 0.35; //作動流体流量(Unit B)[kg/s]



    /*----Constants for Graphics -------------------*/

    var GRAPHDATANUMBER = 60; // Number of Graph Data Points
    var VariableID = 0;
    var PrevioudVariableID = 0;
    var InitialStartValueTemperature = 0, InitialStepWidthTemperature = 4;
    var InitialStartValuePressure = 0, InitialStepWidthPressure = 0.2;
    var InitialStartValuePower = 0, InitialStepWidthPower = 2;
    var haxis=0;

    /* Variables for Graph Data */
    g = new Array(DATANUMBER + 1)
    for(i = 0; i < g.length; i++){
        g[i]= new Array(GRAPHDATANUMBER + 1);
        for(j = 0; j < g[i].length; j++){
            g[i][j] = 0;
        }
    }


    // Conditions for Dynamic Simulation
    var loopNum=0; //Loop Number for Dynamic Simulation
    var DT=0.1;        //Sampling Interval
    var presenttime = 0; // Simulation Time

    // Convergence conditions for Bisection method
    //var Error=0.02;//
    var Error=0.01;// Temperature permissible range. modified by S.Goto 7 March, 2017
    ITERATIONLIMIT=1000; //最適化計算更新回数の制限
    //最適化計算更新回数
    var CNum = 0;
    var ENum = 0;    
    //エネルギーバランス
    var EBe1 = 0.0, EBc1 = 0.0;
    var EBe2 = 0.0, EBc2 = 0.0;
    //蒸発温度,凝縮温度
    var Te1 = 0.0, Tc1 = 0.0;
    var Te2 = 0.0, Tc2 = 0.0;
    //熱交換器の対数平均温度差
    var DTme1 = 0.0, DTmc1 = 0.0;
    var DTme2 = 0.0, DTmc2 = 0.0;

    //タービン出力, ポンプ動力, サイクル出力, サイクル熱効率
    var Wt1 = 0.0, Wt2 = 0.0;
    var Wp1 = 0.0, Wp2 = 0.0;
    var W1 = 0.0, W2 = 0.0;
    var ETA1 = 0.0, ETA2 = 0.0;

    //熱交換量
    var lastQe1 = 0.0, lastQe2 = 0.0;
    var lastQc1 = 0.0, lastQc2 = 0.0;

    /*-------------------------------------*/

    //交換熱量
    var QWS1 = 0.0, QWS2 = 0.0;
    var QCS1 = 0.0, QCS2 = 0.0;
    var Twso1SS = 0.0, Twso2SS = 0.0; //熱交換器出口温度の定常値
    var Tcso1SS = 0.0, Tcso2SS = 0.0; //熱交換器出口温度の定常値
    var QWS1A = 0.0, QWS2A = 0.0; //交換熱量用変数
    var QCS1A = 0.0, QCS2A = 0.0; //交換熱量用変数
    var Twso1 = 0.0, Twso2 = 0.0; //熱交換器出口温度(蒸発器)
    var Tcso1 = 0.0, Tcso2 = 0.0; //熱交換器出口温度(凝縮器)
    var QWS1SS = 0.0, QWS2SS = 0.0; //交換熱量の定常値(蒸発器)
    var QCS1SS = 0.0, QCS2SS = 0.0; //交換熱量の定常値(凝縮器)
    var Twso1A = 0.0, Twso2A = 0.0; //Variables for Warm Seawater Temperature Output
    var Tcso1A = 0.0, Tcso2A = 0.0; //Variables for Warm Seawater Temperature Output
    var MP=0,MPR=0;
    var StepWidth = 0,StartValue = 0;
    var timeoutID = 0, GtimeoutID = 0;


    /***********************************/
    /* Steam Table for Ammonia         */
    /***********************************/

    function STAP(T){ // Saturated Steam Table for Ammonia: Pressure from Temperature 
        var P = 0.0;
        if((0.0 <= T) && (T <= 50.0)){
            P = 4.2950876790616+0.160686057736674*T+0.00233752720738085*Math.pow(T,2)+0.0000157094839995552*Math.pow(T,3)+3.28007873137322*Math.pow(10,-8)*Math.pow(T,4)-3.73500151734348*Math.pow(10,-11)*Math.pow(T,5);
        }else{
            P = NaN;
            alert("T=" + T + ", " + "Pressure: Out of Range.");   // 7 March, 2017 modified by S.Goto
        }
        return P;
    }
    function STALV(T){ // Saturated Steam Table for Ammonia: Liquid Specific Volume from Temperature 
        var V = 0.0;
        if((0.0 <= T) && (T <= 50.0)){
            V = 1.56580224565298+0.00334411338429905*T+0.0000113817018824943*Math.pow(T,2)+1.2021305861484*Math.pow(10,-7)*Math.pow(T,3);
        }else{
            V = NaN;
            alert("T=" + T + ", " + "Volume: Out of Range.");   // 7 March, 2017 modified by S.Goto
        }
        return V;
    }
    function STAVV(T){ // Saturated Steam Table for Ammonia: Vapor Specific Volume from Temperature 
        var V = 0.0;
        if((0.0 <= T) && (T <= 50.0)){
            V = 289.032652116659-10.1306807153039*T+0.199821701739516*Math.pow(T,2)-0.00241395425626903*Math.pow(T,3)+0.0000133107185318683*Math.pow(T,4);
        }else{
            V = NaN;
            alert("T=" + T + ", " + "Volume: Out of Range.");   // 7 March, 2017 modified by S.Goto
        }
        return V;
    }
    function STALH(T){ // Saturated Steam Table for Ammonia: Liquid Specific Enthalpy from Temperature 
        var h = 0.0;
        if((0.0 <= T) && (T <= 50.0)){
            h = -762.749431944304+4.62690087238193*T+0.00334387451429942*Math.pow(T,2)+4.75056071675618*Math.pow(10,-6)*Math.pow(T,3)+1.29349847157388*Math.pow(10,-7)*Math.pow(T,4);
        }else{
            h = NaN;
            alert("STALH: T=" + T + ", " + "Specific Enthalpy: Out of Range.");   // 7 March, 2017 modified by S.Goto
        }
        return h;
    }
    function STAVH(T){ // Saturated Steam Table for Ammonia: Vapor Specific Enthalpy from Temperature 
        var h = 0.0;
        if((0.0 <= T) && (T <= 50.0)){
            h = 499.045936880976+1.06427731926788*T-0.00769053523538193*Math.pow(T,2)-0.0000256653702883702*Math.pow(T,3)-2.34485342639336*Math.pow(10,-7)*Math.pow(T,4);
        }else{
            h = NaN;
            alert("STAVH: T=" + T + ", " + "Specific Enthalpy: Out of Range.");   // 7 March, 2017 modified by S.Goto
        }
        return h;
    }
    function STALS(T){ // Saturated Steam Table for Ammonia: Liquid Specific Entropy from Temperature 
        var s = 0.0;
        if((0.0 <= T) && (T <= 50.0)){
            s = 5.69628720142028+0.0168365424367941*T-0.0000193407494696569*Math.pow(T,2)+3.96028143173146*Math.pow(10,-8)*Math.pow(T,3)+3.48479195859002*Math.pow(10,-10)*Math.pow(T,4);
        }else{
            s = NaN;
            // alert("T=" + T + ", " + "Specific Entropy: Out of Range.");   // 7 March, 2017 modified by S.Goto
        }
        return s;
    }
    function STAVS(T){ // Saturated Steam Table for Ammonia: Vapor Specific Entropy from Temperature 
        var s = 0.0;
        if((0.0 <= T) && (T <= 50.0)){
            s = 10.3156651702603-0.0131082386570281*T+0.0000489161877297273*Math.pow(T,2)-2.91141409528377*Math.pow(10,-7)*Math.pow(T,3);
        }else{
            s = NaN;
            alert("T=" + T + ", " + "Entropy: Out of Range.");   // 7 March, 2017 modified by S.Goto
        }
        return s;
    }
    function STAT(P,s){ // Steam Table for Ammonia: Temperature form Pressure and Specific Enthalpy 
        var a = 0.0, b = 0.0, T = 0.0, SA = 0.0, TA = 0.0, TB = 0.0, TC = 0.0, TD = 0.0;
        SA = s/1000;
        TA = -226.773840+20.2706668*SA+3.43450002*Math.pow(SA,2);//P1.1近似式
        TB = -222.466208264055+18.8385199567052*SA+3.55282683982725*Math.pow(SA,2);//P1.0近似式
        TC = -220.368316547617+18.1276523809518*SA+3.612380952381*Math.pow(SA,2);//P0.9近似式
        TD = -216.361863857194+16.7562835714461*SA+3.72907142856992*Math.pow(SA,2);//P0.8近似式
        //圧力が任意のとき
        if((1000000.0 < P) &&(P <= 1100000.0)){
            a = (1100000.0-P)/100000.0;
            b = 1.0 - a;
            T = TA*b+TB*a;
        }else if((900000.0 < P) &&(P <= 1000000.0)){
            a = (1000000.0-P)/100000.0;
            b = 1.0 - a;
            T = TB*b+TC*a;
        }else if((800000.0 < P) && (P <= 9000000.0)){
            a = (900000.0-P)/100000.0;
            b = 1.0 - a;
            T = TC*b+TD*a;
        }else{
            T = NaN;
            alert("P=" + P + ", " + "s=" + s + ", " + "Temperature: Out of Range.");   // 7 March, 2017 modified by S.Goto
        }
        return T;
    }


    /***********************************/
    /*    計算を行う微分方程式         */
    /*dy/dx= FUNK(xx, yy)              */
    /*戻り値:xx, yyのときのFUNKの値    */
    /***********************************/
    function FUNK(xx, yy, TSS, Tau){
        var yyp;    

        yyp = (TSS-yy)/Tau;   //任意に変更

        return yyp;
    }

    /***********************************************/
    /*       4th order Runge-Kutta method          */
    /*yi: i段目のyの値                             */
    /*xi: i段目のxの値                             */
    /*h : 刻み幅                                   */
    /*戻り値: i+1段目のyの値                       */
    /***********************************************/
    function RungeKutta(yi, xi, h, TSS, Tau){
        var yip = 0.0;
        var w1 = 0.0, w2 = 0.0, w3 = 0.0, w4 = 0.0;

        w1 = h * FUNK(xi,yi,TSS,Tau);
        w2 = h * FUNK(xi+h/2, yi+w1/2,TSS,Tau);
        w3 = h * FUNK(xi+h/2, yi+w2/2,TSS,Tau);
        w4 = h * FUNK(xi+h, yi+w3,TSS,Tau);

        yip = yi + (w1+2*w2+2*w3+w4)/6.0;
        return yip;
    }

    /*--------------ver2 シミュレーション用関数-------------------------------------------------------------------------*/
    function SSimulator2(Twsi1,Twso1,Tcsi1,Tcso1,Twsi2,Twso2,Tcsi2,Tcso2,QWS1,QWS2,QCS1,QCS2){

        var T1,P1,V1,h1,S1,TA,TB,TC,SA,a,b;
        var X,SDD,SD,VDD,VD,HDD,HD;

        //温海水高温熱源
        DTme1 = 1/(UAe1/QWS1);
        //蒸発温度(高温側)
        Te1 = (Math.exp((Twsi1-Twso1)/DTme1) * Twso1 - Twsi1)/(Math.exp((Twsi1 - Twso1)/DTme1) - 1);
        //冷海水高温熱源
        DTmc1 = 1/(UAc1/QCS1);
        // alert("Ver.2.0: "+ DTmc1+UAc1+QCS1);
        //凝縮温度(高温側)
        Tc1 = (Tcso1 - Math.exp((Tcsi1-Tcso1)/DTmc1) * Tcsi1)/(1 - Math.exp((Tcsi1 - Tcso1)/DTmc1));
        //温海水低温熱源
        DTme2 = 1/(UAe2/QWS2);
        //蒸発温度(低温側)
        Te2 = (Math.exp((Twsi2-Twso2)/DTme2) * Twso2 - Twsi2)/(Math.exp((Twsi2 - Twso2)/DTme2) - 1);
        //冷海水低温熱源
        DTmc2 = 1/(UAc2/QCS2);
        //凝縮温度(低温側)
        Tc2 = (Tcso2 - Math.exp((Tcsi2-Tcso2)/DTmc2) * Tcsi2)/(1 - Math.exp((Tcsi2 - Tcso2)/DTmc2));

        //------------------------------//
        //蒸発器部(Point 1, 4B)
        T1 = Te1;  //Point 1の温度は蒸発温度と同じ
        P1 = STAP(T1);
        P1 = P1*Math.pow(10,5);
        V1 = STAVV(T1);
        V1 = V1*Math.pow(10,-3);
        h1 = STAVH(T1);
        h1 = h1*1000;  //単位変換J/g-->J/kg
        S1 = STAVS(T1);
        S1 = S1*1000;  //単位変換J/g-->J/kg
        T[0]=T1;
        P[0]=P1;
        V[0]=V1;
        h[0]=h1;
        S[0]=S1;
        //-------------------------------//
        //凝縮器部(Point 3)
        T1 = Tc1;   //Point 3の温度は凝縮温度と同じ
        P1 = STAP(T1);
        P1 = P1*Math.pow(10,5);  //単位変換 bar-->Pa
        V1 = STALV(T1);
        V1 = V1*Math.pow(10,-3);   //単位変換
        // alert("Point3: "+T1);
        h1 = STALH(T1);
        h1 = h1*1000;  //単位変換
        S1 = STALS(T1);
        S1 = S1*1000;  //単位変換
        T[2] = T1;
        P[2] = P1;
        V[2] = V1;
        h[2] = h1;
        S[2] = S1;
        //-------------------------------//
        //Point 2(湿り蒸気)
        S1 = S[0];   //Point 1->2は等エントロピー変化
        T1 = Tc1;     //凝縮温度と同じ
        P1 = P[2];   //Point 2, Point 3は等圧力
        //乾き度の計算から修正20141103
        SDD = STAVS(T1);
        SDD = SDD*1000;  //単位変換J/g-->J/kg
        SD = STALS(T1);
        SD = SD*1000;  //単位変換J/g-->J/kg
        X = (S1-SD)/(SDD-SD);
        // alert("Point2: "+T1);
        HDD = STAVH(T1);
        HD = STALH(T1);
        h1 = (HD+X*(HDD-HD));
        h1 = h1*1000;  //単位変換
        VDD = STAVV(T1);
        VD = STALV(T1);
        V1 = (VD+X*(VDD-VD));
        V1 = V1*Math.pow(10,-3);   //単位変換
        T[1] = T1;
        P[1] = P1;
        V[1] = V1;
        h[1] = h1;
        S[1] = S1;
        //-------------------------------//
        //Point 4
        P1 = P[0];    //Point 4とPoint 1は等圧
        S1 = S[2];    //Point 4とPoint 3は等エントロピー
        T1 = STAT(P1,S1);
        V1 = V[2];    //Point 4とPoint 3の比体積はほとんど変わらないため
        h1 = h[2]+V1*(P1-P[2]);   //Point 4の比エンタルピーの計算

        T[3] = T1;
        P[3] = P1;
        V[3] = V1;
        h[3] = h1;
        S[3] = S1;

        //---------------------------------------//

        //蒸発器部(Point 5, 8B)
        T1 = Te2;  //Point 5の温度は蒸発温度と同じ
        P1 = STAP(T1);
        P1 = P1*Math.pow(10,5);
        V1 = STAVV(T1);
        V1 = V1*Math.pow(10,-3);
        h1 = STAVH(T1);
        h1 = h1*1000;  //単位変換J/g-->J/kg
        S1 = STAVS(T1);
        S1 = S1*1000;  //単位変換J/g-->J/kg
        T[4]=T1;
        P[4]=P1;
        V[4]=V1;
        h[4]=h1;
        S[4]=S1;

        //---------------------------------------//

        //凝縮器部(Point 7)
        T1 = Tc2;   //Point 7の温度は凝縮温度と同じ
        P1 = STAP(T1);
        P1 = P1*Math.pow(10,5);  //単位変換 bar-->Pa
        V1 = STALV(T1);
        V1 = V1*Math.pow(10,-3);   //単位変換
        // alert("Point7: "+T1);
        h1 = STALH(T1);
        h1 = h1*1000;  //単位変換
        S1 = STALS(T1);
        S1 = S1*1000;  //単位変換
        T[6] = T1;
        P[6] = P1;
        V[6] = V1;
        h[6] = h1;
        S[6] = S1;

        //---------------------------------------//

        //Point 6(湿り蒸気)
        S1 = S[4];   //Point 5->6は等エントロピー変化
        T1 = Tc2;     //凝縮温度と同じ
        P1 = P[6];   //Point 6, Point 7は等圧力
        //乾き度の計算から修正20141103
        SDD = STAVS(T1);
        SDD = SDD*1000;  //単位変換J/g-->J/kg
        SD = STALS(T1);
        SD = SD*1000;  //単位変換J/g-->J/kg
        X = (S1-SD)/(SDD-SD);
        // alert("Point6: "+T1);
        HDD = STAVH(T1);
        HD = STALH(T1);
        h1 = (HD+X*(HDD-HD));
        h1 = h1*1000;  //単位変換
        VDD = STAVV(T1);
        VD = STALV(T1);
        V1 = (VD+X*(VDD-VD));
        V1 = V1*Math.pow(10,-3);   //単位変換
        T[5] = T1;
        P[5] = P1;
        V[5] = V1;
        h[5] = h1;
        S[5] = S1;
        //---------------------------------------//

        //Point 8
        P1 = P[4];    //Point 8とPoint 5は等圧
        S1 = S[6];    //Point 8とPoint 7は等エントロピー
        T1 = STAT(P1,S1);
        V1 = V[6];    //Point 8とPoint 7の比体積はほとんど変わらないため
        h1 = h[6]+V1*(P1-P[6]);   //Point 8の比エンタルピーの計算

        T[7] = T1;
        P[7] = P1;
        V[7] = V1;
        h[7] = h1;
        S[7] = S1;
        //---------------------------------------//
        //単位変換
        ////PをMPaへ
        P[0] = P[0]/1000000;
        P[1] = P[1]/1000000;
        P[2] = P[2]/1000000;
        P[3] = P[3]/1000000;
        P[4] = P[4]/1000000;
        P[5] = P[5]/1000000;
        P[6] = P[6]/1000000;
        P[7] = P[7]/1000000;
        ////JをKJへ
        h[0] = h[0]/1000;
        h[1] = h[1]/1000;
        h[2] = h[2]/1000;
        h[3] = h[3]/1000;
        h[4] = h[4]/1000;
        h[5] = h[5]/1000;
        h[6] = h[6]/1000;
        h[7] = h[7]/1000;
        ////J/KgkをKJ/KgKへ
        S[0] = S[0]/1000;
        S[1] = S[1]/1000;
        S[2] = S[2]/1000;
        S[3] = S[3]/1000;
        S[4] = S[4]/1000;
        S[5] = S[5]/1000;
        S[6] = S[6]/1000;
        S[7] = S[7]/1000;
        //---------------------------------------//
    }

    /*--------------------Ver2.5------------------------------*/
    function SSimulator25(Twsi1, Mws, Tcsi2, Mcs, Mwf1, Mwf2){
        var ITwso2=0,ITcso1=0,ITwso1=0,ITcso2=0;
        var ERR1=0,ERR2=0,ERR3=0,ERR4=0;
        var bc1=0.0,bc2=0.0,aw1=0.0,aw2=0.0;
        var ac1=0.0,ac2=0.0,abc1=0.0,abc2=0.0,aw1=0.0,aw2=0.0;
        var bw1=0.0;
        var h1a=0.0,h2a=0.0,h0a=0.0,h3a=0.0,h4a=0.0,h6a=0.0,h5a=0.0,h7a=0.0;
        var T1=0, P1=0, V1=0, h1=0, S1=0,TA=0,TB=0,TC=0,SA=0,a=0,b=0;
        var X=0,SDD=0,SD=0,VDD=0,VD=0,HDD=0,HD=0;
        var IniTempDiff=5.0; //Very Very Important !! The convergence conditions highly depend on the IniTempDiff parameter. 収束の範囲に大きく影響を与える． 14 March, 2017
        var EC=0.0,EC1=0.0,EC2=0.0,EE=0.0,EE1=0.0,EE2=0.0,ECC=0.0,ECC1=0.0,ECC2=0.0,EEE=0.0,EEE1=0.0,EEE2=0.0; 
        //-------------------------------//
        /*Tcso2(冷海水低温側出口)の仮定*/
        ITcso2=0;
        bc2 = Tcsi2+IniTempDiff;
        while(ITcso2<ITERATIONLIMIT){
            ITcso2++; 
            // alert("Cp: "+Cp+bc2+Tcsi2);
            QCS2SS = Mcs*Cp*(bc2-Tcsi2);   //冷海水低温熱源の熱交換量
            DTmc2 = 1/(UAc2/QCS2SS);
            //凝縮温度(低温側)
            Tc2 = (bc2 - Math.exp((Tcsi2-bc2)/DTmc2) * Tcsi2)/(1 - Math.exp((Tcsi2-bc2)/DTmc2));

            //-------------------------------//
            /*Twso1(温海水高温側出口)の仮定*/
            ITwso1=0;
            aw1 = Twsi1-IniTempDiff;
            while(ITwso1<ITERATIONLIMIT){
                ITwso1++; 
                QWS1SS = Mws*Cp*(Twsi1-aw1);   //温海水低温熱源の熱交換量
                DTme1 = 1/(UAe1/QWS1SS);
                //蒸発温度(高温側)
                Te1 = (Math.exp((Twsi1-aw1)/DTme1) * aw1 - Twsi1)/(Math.exp((Twsi1-aw1)/DTme1) - 1);

                //------------------------------//
                /*Tcso1(冷海水高温側出口)の仮定*/
                ITcso1=0;
                ac1 = Tcsi1;
                bc1 = bc2+IniTempDiff;
                while(ITcso1<ITERATIONLIMIT){
                    ITcso1++; 
                    QCS1SS = Mcs*Cp*(bc1-bc2);   //冷海水高温熱源の熱交換量
                    DTmc1 = 1/(UAc1/QCS1SS);
                    //凝縮温度(高温側)
                    Tc1 = (bc1 - Math.exp((bc2-bc1)/DTmc1) * bc2)/(1 - Math.exp((bc2 - bc1)/DTmc1));
                    // alert("Ver.2.5: "+ Tc1+bc1+bc2+DTmc1+UAc1+QCS1SS+Mcs+Cp);
                    //------------------------------//
                    //蒸発器部(Point 1, 4B)
                    T1 = Te1;  //Point 1の温度は蒸発温度と同じ
                    //****近似式
                    P1 = STAP(T1);
                    P1 = P1*Math.pow(10,5);
                    V1 = STAVV(T1);
                    V1 = V1*Math.pow(10,-3);
                    h1 = STAVH(T1);
                    h1 = h1*1000;  //単位変換J/g-->J/kg
                    S1 = STAVS(T1);
                    S1 = S1*1000;  //単位変換J/g-->J/kg
                    T[0]=T1;
                    P[0]=P1;
                    V[0]=V1;
                    h[0]=h1;
                    S[0]=S1;
                    //-------------------------------//
                    //凝縮器部(Point 3)
                    T1 = Tc1;   //Point 3の温度は凝縮温度と同じ
                    P1 = STAP(T1);
                    P1 = P1*Math.pow(10,5);  //単位変換 bar-->Pa
                    V1 = STALV(T1);
                    V1 = V1*Math.pow(10,-3);   //単位変換
                    // alert("Ver.2.5 Point3: "+ T1);
                    h1 = STALH(T1);
                    h1 = h1*1000;  //単位変換
                    S1 = STALS(T1);
                    S1 = S1*1000;  //単位変換
                    T[2] = T1;
                    P[2] = P1;
                    V[2] = V1;
                    h[2] = h1;
                    S[2] = S1;
                    //-------------------------------//
                    //Point 2(湿り蒸気)
                    S1 = S[0];   //Point 1->2は等エントロピー変化
                    T1 = Tc1;     //凝縮温度と同じ
                    P1 = P[2];   //Point 2, Point 3は等圧力
                    //乾き度の計算から修正20141103
                    SDD = STAVS(T1);
                    SDD = SDD*1000;  //単位変換J/g-->J/kg
                    SD = STALS(T1);
                    SD = SD*1000;  //単位変換J/g-->J/kg
                    X = (S1-SD)/(SDD-SD);
                    // alert("Ver.2.5 Point2: "+T1);
                    HDD = STAVH(T1);
                    HD = STALH(T1);
                    h1 = (HD+X*(HDD-HD));
                    h1 = h1*1000;  //単位変換
                    VDD = STAVV(T1);
                    VD = STALV(T1);
                    V1 = (VD+X*(VDD-VD));
                    V1 = V1*Math.pow(10,-3);   //単位変換
                    T[1] = T1;
                    P[1] = P1;
                    V[1] = V1;
                    h[1] = h1;
                    S[1] = S1;
                    //-------------------------------//
                    //Point 4
                    P1 = P[0];    //Point 4とPoint 1は等圧
                    S1 = S[2];    //Point 4とPoint 3は等エントロピー
                    T1 = STAT(P1,S1);
                    V1 = V[2];    //Point 4とPoint 3の比体積はほとんど変わらないため
                    h1 = h[2]+V1*(P1-P[2]);   //Point 4の比エンタルピーの計算
                    T[3] = T1;
                    P[3] = P1;
                    V[3] = V1;
                    h[3] = h1;
                    S[3] = S1;
                    //-------------------------------//

                    // エネルギーバランスの値の確認
                    //****Tcso1の更新****

                    EC1=1.0-((Mcs*Cp*(Tcso1SS-bc2))/(Mwf1*(h[1]-h[2])));
                    EC2=1.0-((Mcs*Cp*(ac1-bc2))/(Mwf1*(h1a-h2a)));
                    //  alert("EC1=" + EC1 + " EC2=" + EC2 +" Tcso1SS=" + Tcso1SS +" Mcs="+Mcs+" Mwf1="+Mwf1+" h[1]="+h[1]+" h[2]="+h[2]+" Iteration="+ ITcso1);   //
                    if(ITcso1==1){  //1回目
                        h1a = h[1];
                        h2a = h[2];
                        Tcso1SS=bc2;
                        ac1=Tcso1SS;
                    }else{
                        if(EC1*EC2<0){
                            //  alert("EC1*EC2<0: EC1=" + EC1+" EC2=" + EC2+" Tcso1SS=" + Tcso1SS +" Iteration="+ ITcso1+" ac1="+ac1+" bc1="+bc1);   //
                            bc1 = Tcso1SS;
                        }else{
                            //  alert("EC1*EC2>=0: EC1=" + EC1+" EC2=" + EC2+" Tcso1SS=" + Tcso1SS +" Iteration="+ ITcso1+" ac1="+ac1+" bc1="+bc1);   //
                            ac1 = Tcso1SS;  //ac1を更新
                            h1a = h[1]; //ac1でのエンタルピーを更新
                            h2a = h[2];
                        }
                        if(Math.abs(EC1)<Error){
                            ERR1 = EC1;   //エネルギーバランスの計算
                            break;
                        }
                    }
                    Tcso1SS = (ac1+bc1)/2.0;   //Tcsoの更新
                    //     alert("ITcso1="+ITcso1+" Tcso1SS=" + Tcso1SS + " ERR1=" + ERR1);
                    //****Tcso1の更新ここまで****
                } //ここまでの範囲で凝縮器のエネルギーバランスを最適化
                //  alert("ERR1=" + EC1+" Tcso1SS=" + Tcso1SS +" Iteration="+ ITcso1);   //
                bc1=Tcso1SS;
                CNum = ITcso1; 
                EBc1=ERR1;
                ERR1=0.0;
                //----------------------------------------//
                //****Twso1の更新****
                //****二分法(20141012追加)

                EE1=1.0-((Mws*Cp*(Twsi1-Twso1SS))/(Mwf1*(h[0]-h[3])));
                EE2=1.0-((Mws*Cp*(Twsi1-aw1))/(Mwf1*(h0a-h3a)));
                //  alert("EE1=" + EE1 + " EE2=" + EE2 +" Twso1SS=" + Twso1SS +" Mws="+Mws+" Mwf1="+Mwf1+" h[0]="+h[0]+" h[3]="+h[3]+" Iteration="+ ITwso1);   //
                if(ITwso1==1){  //1回目
                    h0a = h[0];
                    h3a = h[3];
                    Twso1SS=Twsi1;
                    bw1 = Twso1SS;
                    //  alert(" Iteration="+ ITwso1+" EE1=" + EE1+" EE2=" + EE2+" Twso1SS=" + Twso1SS+" aw1="+aw1+" bw1=" + bw1);   //
                }else{
                    if(EE1*EE2<0){
                        //  alert("EE1*EE2<0: EE1=" + EE1+" EE2=" + EE2+" Twso1SS=" + Twso1SS +" Iteration="+ ITwso1+" aw1="+aw1+" bw1="+bw1);   //
                        bw1 = Twso1SS;
                    }else{
                        //  alert("EE1*EE2>=0: EE1=" + EE1+" EE2=" + EE2+" Twso1SS=" + Twso1SS +" Iteration="+ ITwso1+" aw1="+aw1+" bw1="+bw1);   //
                        aw1 = Twso1SS;  //ac1を更新
                        h0a = h[0]; //ac1でのエンタルピーを更新
                        h3a = h[3];
                    }
                    if(Math.abs(EE1)<Error){
                        ERR2 = EE1;   //エネルギーバランスの計算
                        //  alert("ERR2=" + EE1+" Twso1SS=" + Twso1SS +" Iteration="+ ITwso1);   //
                        break;
                    }
                }
                if(ITwso1>=ITERATIONLIMIT){
                    alert("EE1=" + EE1 + " EE2=" + EE2 +  " aw1=" + aw1 +  " ITwso1="+ ITwso1 + " Twso1 Not converged.");   //1000回実行して収束しないとき異常終了 7 March, 2017 modified by S.Goto 
                }
                Twso1SS = (aw1+bw1)/2.0;   //Tcsoの更新
                //     alert("ITwso1="+ITwso1+" Twso1SS=" + Twso1SS +" aw1="+aw1+" bw1="+bw1+ " ERR2=" + ERR2);
                //****Twso1の更新ここまで****
            }//ここまでの範囲で凝縮器のエネルギーバランスを最適化
            //  alert("Twso1SS=" + Twso1SS +" Iteration="+ ITwso1);   //
            aw1=Twso1SS;
            CNum = ITwso1; 
            EBe1=ERR2;
            ERR2=0.0;
            //---------------------------------------//
            /*Twso2(温海水低温側出口)の仮定*/
            ITwso2=0;
            aw2 =aw1-IniTempDiff;
            while(ITwso2<ITERATIONLIMIT){
                ITwso2++; 
                //****Ver2.5追加箇所
                QWS2SS = Mws*Cp*(aw1-aw2);   //温海水低温熱源の熱交換量
                DTme2 = 1/(UAe2/QWS2SS);
                //****Ver2.5追加箇所ここまで
                //蒸発温度(低温側)
                Te2 = (Math.exp((aw1-aw2)/DTme2) * aw2 - aw1)/(Math.exp((aw1-aw2)/DTme2) - 1);
                //---------------------------------------//
                //蒸発器部(Point 5, 8B)
                T1 = Te2;  //Point 5の温度は蒸発温度と同じ
                P1 = STAP(T1);
                P1 = P1*Math.pow(10,5);
                V1 = STAVV(T1);
                V1 = V1*Math.pow(10,-3);
                h1 = STAVH(T1);
                h1 = h1*1000;  //単位変換J/g-->J/kg
                S1 = STAVS(T1);
                S1 = S1*1000;  //単位変換J/g-->J/kg
                T[4]=T1;
                P[4]=P1;
                V[4]=V1;
                h[4]=h1;
                S[4]=S1;
                //---------------------------------------//
                //凝縮器部(Point 7)
                T1 = Tc2;   //Point 7の温度は凝縮温度と同じ
                P1 = STAP(T1);
                P1 = P1*Math.pow(10,5);  //単位変換 bar-->Pa
                V1 = STALV(T1);
                V1 = V1*Math.pow(10,-3);   //単位変換
                // alert("Ver. 2.5 Point7: "+T1);
                h1 = STALH(T1);
                h1 = h1*1000;  //単位変換
                S1 = STALS(T1);
                S1 = S1*1000;  //単位変換     
                T[6] = T1;
                P[6] = P1;
                V[6] = V1;
                h[6] = h1;
                S[6] = S1;
                //---------------------------------------//
                //Point 6(湿り蒸気)
                S1 = S[4];   //Point 5->6は等エントロピー変化
                T1 = Tc2;     //凝縮温度と同じ
                P1 = P[6];   //Point 6, Point 7は等圧力
                //乾き度の計算から修正20141103
                SDD = STAVS(T1);
                SDD = SDD*1000;  //単位変換J/g-->J/kg
                SD = STALS(T1);
                SD = SD*1000;  //単位変換J/g-->J/kg
                X = (S1-SD)/(SDD-SD);
                // alert("Ver. 2.5 Point6: "+T1);
                HDD = STAVH(T1);
                HD = STALH(T1);
                h1 = (HD+X*(HDD-HD));
                h1 = h1*1000;  //単位変換
                VDD = STAVV(T1);
                VD = STALV(T1);
                V1 = (VD+X*(VDD-VD));
                V1 = V1*Math.pow(10,-3);   //単位変換
                T[5] = T1;
                P[5] = P1;
                V[5] = V1;
                h[5] = h1;
                S[5] = S1;
                //---------------------------------------//
                //Point 8
                P1 = P[4];    //Point 8とPoint 5は等圧
                S1 = S[6];    //Point 8とPoint 7は等エントロピー
                T1 = STAT(P1,S1);
                V1 = V[6];    //Point 8とPoint 7の比体積はほとんど変わらないため
                h1 = h[6]+V1*(P1-P[6]);   //Point  8の比エンタルピーの計算
                T[7] = T1;
                P[7] = P1;
                V[7] = V1;
                h[7] = h1;
                S[7] = S1;
                //---------------------------------------//
                //****Twso2の更新
                //****二分法(20141012追加)
                EEE1=1.0-((Mws*Cp*(aw1-Twso2SS))/(Mwf2*(h[4]-h[7])));
                EEE2=1.0-((Mws*Cp*(aw1-aw2))/(Mwf2*(h4a-h7a)));
                //  alert("EEE1=" + EEE1 + " EEE2=" + EEE2 +" Twso2SS=" + Twso2SS +" Mws="+Mws+" Mwf2="+Mwf2+" h[4]="+h[4]+" h[7]="+h[7]+" Iteration="+ ITwso2);   //
                if(ITwso2==1){  //1回目
                    h4a = h[4];
                    h7a = h[7];
                    Twso2SS=aw1;
                    bw2 = Twso2SS;
                    //  alert(" Iteration="+ ITwso2+" EEE1=" + EEE1+" EEE2=" + EEE2+" Twso2SS=" + Twso2SS+" aw2="+aw2+" bw2=" + bw2);   //
                    }else{
                        if(EEE1*EEE2<0){
                            //  alert("EEE1*EEE2<0: EEE1=" + EEE1+" EEE2=" + EEE2+" Twso2SS=" + Twso2SS +" Iteration="+ ITwso2+" aw2="+aw2+" bw2="+bw2);   //
                            bw2 = Twso2SS;
                        }else{
                            //  alert("EEE1*EEE2>=0: EEE1=" + EEE1+" EEE2=" + EEE2+" Twso2SS=" + Twso2SS +" Iteration="+ ITwso2+" aw2="+aw2+" bw2="+bw2);   //
                            aw2 = Twso2SS;  //ac1を更新
                            h4a = h[4]; //aw2でのエンタルピーを更新
                            h7a = h[7];
                        }
                        if(Math.abs(EEE1)<Error){
                            ERR3 = EEE1;   //エネルギーバランスの計算
                            //  alert("ERR3=" + EEE1+" Twso2SS=" + Twso2SS +" Iteration="+ ITwso2);   //
                            break;
                        }
                    }
                    if(ITwso2>=ITERATIONLIMIT){
                        alert("EEE1=" + EEE1 + " EEE2=" + EEE2 +  " aw2=" + aw2 +  " ITwso2="+ ITwso2 + " Twso2 Not converged.");   //1000回実行して収束しないとき異常終了 7 March, 2017 modified by S.Goto 
                    }
                    Twso2SS = (aw2+bw2)/2.0;   //Tcsoの更新
                    //     alert("ITwso2="+ITwso2+" Twso2SS=" + Twso2SS +" aw2="+aw2+" bw2="+bw2+ " ERR3=" + ERR3);
                    //****Twso2の更新ここまで****
                }//ここまでの範囲で蒸発器のエネルギーバランスを最適化
                //  alert("Twso2SS=" + Twso2SS +" Iteration="+ ITwso2);   //
                aw2=Twso2SS;
                CNum = ITwso2; 
                EBe2=ERR3;
                ERR3=0.0;

                //---------------------------------------//
                //****Tcso2の更新
                //****二分法(20141012追加)
                ECC1=1.0-((Mcs*Cp*(Tcso2SS-Tcsi2))/(Mwf2*(h[5]-h[6])));
                ECC2=1.0-((Mcs*Cp*(ac2-Tcsi2))/(Mwf2*(h5a-h6a)));
                //  alert("ECC1=" + ECC1 + " ECC2=" + ECC2 +" Tcso2SS=" + Tcso2SS +" Mcs="+Mcs+" Mwf2="+Mwf2+" h[5]="+h[5]+" h[6]="+h[6]+" Iteration="+ ITcso2);   //
                if(ITcso2==1){  //1回目
                    h5a = h[5];
                    h6a = h[6];
                    Tcso2SS=Tcsi2;
                    ac2=Tcso2SS;
                }else{
                    if(ECC1*ECC2<0){
                        //  alert("ECC1*ECC2<0: ECC1=" + ECC1+" ECC2=" + ECC2+" Tcso2SS=" + Tcso2SS +" Iteration="+ ITcso2+" ac2="+ac2+" bc2="+ bc2);   //
                        bc2 = Tcso2SS;
                    }else{
                        //  alert("ECC1*ECC2>=0: ECC1=" + ECC1+" ECC2=" + ECC2+" Tcso2SS=" + Tcso2SS +" Iteration="+ ITcso2+" ac2="+ac2+" bc2="+ bc2);   //
                        ac2 = Tcso2SS;  //ac2を更新
                        h5a = h[5]; //ac2でのエンタルピーを更新
                        h6a = h[6];
                    }
                    if(Math.abs(ECC1)<Error){
                        ERR4 = ECC1;   //エネルギーバランスの計算
                        break;
                    }
                }
                Tcso2SS = (ac2+bc2)/2.0;   //Tcsoの更新
                //     alert("ITcso2="+ITcso2+" Tcso2SS=" + Tcso2SS + " ERR4=" + ERR4);
                //****Tcso2の更新ここまで****
            } //ここまでの範囲で凝縮器のエネルギーバランスを最適化
            //  alert("ERR2=" + ECC1+" Tcso2SS=" + Tcso2SS +" Iteration="+ ITcso2);   //
            bc2=Tcso2SS;
            CNum = ITcso2; 
            EBc2=ERR4;
            ERR4=0.0;
            //---------------------------------------//
            //単位変換
            ////PをMPaへ
            P[0] = P[0]/1000000;
            P[1] = P[1]/1000000;
            P[2] = P[2]/1000000;
            P[3] = P[3]/1000000;
            P[4] = P[4]/1000000;
            P[5] = P[5]/1000000;
            P[6] = P[6]/1000000;
            P[7] = P[7]/1000000;

            ////JをKJへ
            h[0] = h[0]/1000;
            h[1] = h[1]/1000;
            h[2] = h[2]/1000;
            h[3] = h[3]/1000;
            h[4] = h[4]/1000;
            h[5] = h[5]/1000;
            h[6] = h[6]/1000;
            h[7] = h[7]/1000;

            ////J/KgkをKJ/KgKへ
            S[0] = S[0]/1000;
            S[1] = S[1]/1000;
            S[2] = S[2]/1000;
            S[3] = S[3]/1000;
            S[4] = S[4]/1000;
            S[5] = S[5]/1000;
            S[6] = S[6]/1000;
            S[7] = S[7]/1000;
    }



    //入力値の読み込み
    function InputData(){
        var tw=document.getElementById("tw");
        var tc=document.getElementById("tc");
        var mw=document.getElementById("mw");
        var mc=document.getElementById("mc");
        var mwf1=document.getElementById("mwf1");
        var mmf2=document.getElementById("mwf2");

        Twsi1 = parseFloat(tw.value);
        Tcsi2 = parseFloat(tc.value);
        Mws = parseFloat(mw.value);
        Mcs = parseFloat(mc.value);
        Mwf1 = parseFloat(mwf1.value);
        Mwf2 = parseFloat(mwf2.value);

        // alert("Input Completed")
        Sleep( 500 ); // 500ms sleeep
        document.getElementById("tw").style.color="blue";
        document.getElementById("tc").style.color="blue";
        document.getElementById("mw").style.color="blue";
        document.getElementById("mc").style.color="blue";
        document.getElementById("mwf1").style.color="blue";
        document.getElementById("mwf2").style.color="blue";
    }

    function Sleep( T ){ 
        var d1 = new Date().getTime(); 
        var d2 = new Date().getTime(); 
        while( d2 < d1+T ){    //T秒待つ 
            d2=new Date().getTime(); 
        } 
        return; 
    } 

    //******動的シミュレータ
    function DSimulator(){

        //初期化
        Twso1SS = 0; 
        Tcso1SS = 0;
        Twso2SS = 0;
        Tcso2SS = 0;

        SSimulator25(Twsi1, Mws, Tcsi2, Mcs, Mwf1, Mwf2);    //静的シミュレータ Ver2.5

        // Time constants of evaporator and condenser
        var TauE1 = cw*rhow*Ae1*sigmae1/(cw*Mws+Ue1*Ae1);//温海水出口温度の時定数 (UnitA) [s]
        var TauE2 = cw*rhow*Ae2*sigmae2/(cw*Mws+Ue2*Ae2);//温海水出口温度の時定数 (UnitB) [s]
        var TauC1 = cw*rhow*Ac1*sigmac1/(cw*Mcs+Uc1*Ac1);//温海水出口温度の時定数 (UnitA) [s]
        var TauC2 = cw*rhow*Ac2*sigmac2/(cw*Mcs+Uc2*Ac2);//温海水出口温度の時定数 (UnitB) [s]
        var TauEQ1 = ca*rhoa*Ae1*sigmae1/(ca*Mwf1+Ue1*Ae1);//温海水出口温度の時定数 (UnitA) [s]
        var TauEQ2 = ca*rhoa*Ae2*sigmae2/(ca*Mwf2+Ue2*Ae2);//温海水出口温度の時定数 (UnitB) [s]
        var TauCQ1 = ca*rhoa*Ac1*sigmac1/(ca*Mwf1+Uc1*Ac1);//温海水出口温度の時定数 (UnitA) [s]
        var TauCQ2 = ca*rhoa*Ac2*sigmac2/(ca*Mwf2+Uc2*Ac2);//温海水出口温度の時定数 (UnitB) [s]

        // alert("Time Constants;\n" +  "TauE1=" + eval(TauE1) + " TauE2=" + eval(TauE2) + " TauC1=" + eval(TauC1) + " TauC2=" + eval(TauC2) + "\n" + "TauEQ1=" + eval(TauEQ1) + " TauEQ2=" + eval(TauEQ2) + " TauCQ1=" + eval(TauCQ1) + " TauCQ2=" + eval(TauCQ2) + "\n"); 


        /*--------------熱交換器出口温度の計算(初期値から定常値まで1次系で変化すると仮定)--------------*/
        if(loopNum==0){
            Twso1 = Twso1SS;
            Twso2 = Twso2SS;
            Tcso1 = Tcso1SS;
            Tcso2 = Tcso2SS;
        }else{
            Twso1 = RungeKutta(Twso1A, loopNum*DT, DT, Twso1SS, TauE1); //ルンゲクッタ(蒸発器1出口温度 UnitA)    
            Twso2 = RungeKutta(Twso2A, loopNum*DT, DT, Twso2SS, TauE2); //ルンゲクッタ(蒸発器2出口温度 UnitB)
            Tcso1 = RungeKutta(Tcso1A, loopNum*DT, DT, Tcso1SS, TauC1); //ルンゲクッタ(凝縮器1出口温度 UnitA)
            Tcso2 = RungeKutta(Tcso2A, loopNum*DT, DT, Tcso2SS, TauC2); //ルンゲクッタ(凝縮器2出口温度 UnitB)
        }
        Twso1A = Twso1;
        Twso2A = Twso2;
        Tcso1A = Tcso1; 
        Tcso2A = Tcso2;   //次段でルンゲ・クッタに使用するデータ

        //****熱交換器出口温度の計算ここまで
        //****交換熱量の計算（交換熱量の変化が1次系と仮定）20141203
        if(loopNum==0){
            QWS1 = QWS1SS;   //高温熱源の熱交換量
            QWS2 = QWS2SS;   //高温熱源の熱交換量
            QCS1 = QCS1SS;   //冷海水が作動流体から受け取る熱量
            QCS2 = QCS2SS;   //冷海水が作動流体から受け取る熱量
        }else{
            QWS1 = RungeKutta(QWS1A, loopNum*DT, DT, QWS1SS, TauEQ1); //ルンゲ・クッタ(蒸発器の交換熱量 UnitA)
            QWS2 = RungeKutta(QWS2A, loopNum*DT, DT, QWS2SS, TauEQ2); //ルンゲ・クッタ(蒸発器の交換熱量 UnitB)
            QCS1 = RungeKutta(QCS1A, loopNum*DT, DT, QCS1SS, TauCQ1); //ルンゲ・クッタ(凝縮器の交換熱量 UnitA)
            QCS2 = RungeKutta(QCS2A, loopNum*DT, DT, QCS2SS, TauCQ2); //ルンゲ・クッタ(凝縮器の交換熱量 UnitB)
        }
        QWS1A=QWS1;
        QWS2A=QWS2;
        QCS1A=QCS1;
        QCS2A=QCS2;  //次段でルンゲ・クッタに使用するデータ
        //****交換熱量の計算ここまで

        Twsi2=Twso1;
        Tcsi1=Tcso2;

        SSimulator2(Twsi1,Twso1,Tcsi1,Tcso1,Twsi2,Twso2,Tcsi2,Tcso2,QWS1,QWS2,QCS1,QCS2); //Static Simulation Ver.2

        /*---------------------------------------------------------------------------------------------*/


        ////タービン出力(単位はkW)
        Wt1 = etaT*Mwf1*(h[0]-h[1]);
        Wt2 = etaT*Mwf2*(h[4]-h[5]);

        ////作動流体ポンプ動力(単位はkW)
        Wp1 = Mwf1*(h[3]-h[2])/etaP;
        Wp2 = Mwf2*(h[7]-h[6])/etaP;
        ////サイクルの出力(単位はkW)
        W1 = etaG*Wt1-Wp1;
        W2 = etaG*Wt2-Wp2;
        ///蒸発器内交換熱量Qe
        lastQe1 = Mwf1*(h[0]-h[3]);
        lastQe2 = Mwf2*(h[4]-h[7]);
        ///凝縮器内交換熱量Qc
        lastQc1 = Mwf1*(h[1]-h[2]);
        lastQc2 = Mwf2*(h[5]-h[6]);
        ////サイクルの熱交換効率(単位は%)
        ETA1 = 100*W1/lastQe1;
        ETA2 = 100*W2/lastQe2;

        /* 最大仕事量 */
        MP = (Mws*Cp/1000*Mcs*Cp/1000)/Math.pow(Number(Mws*Cp/1000+Mcs*Cp/1000),2)*( (2*Mws*Cp/1000+Mcs*Cp/1000)*(Number(Twsi1)+273)-3*(Mws*Cp/1000*Math.pow( (Number(Twsi1)+273),2/3)*Math.pow( (Number(Tcsi2)+273),1/3)+Mcs*Cp/1000*Math.pow( (Number(Twsi1)+273),1/3)*Math.pow( (Number(Tcsi2)+273),2/3))+(Mws*Cp/1000+2*Mcs*Cp/1000)*(Number(Tcsi2)+273) );
        /* 最大仕事率 */
        MPR = ( ( Mws*Cp/1000*( (Number(Twsi1)+273)-( Number(Twso2)+273))-Mcs*Cp/1000*( (Number(Tcso1)+273)-(Number(Tcsi2)+273) ) )/ MP )*100;

        loopNum++;   //1段進める

        return;
     }
     //*****動的シミュレータここまで


    //Main Calculation for Simulator//
    function mt(){

        DSimulator();   //動的シミュレータ
     
        if(loopNum%10==0){

            // Data Output of Simulation Results for Graphic and Tables in Simulator
            //温海水入口温度
            document.getElementById("data1").innerHTML = Twsi1.toFixed(2);
            tab1.rows[0].cells[1].innerText = (Number(data1.innerHTML)).toFixed(2);
            //冷海水入口温度
            document.getElementById("data2").innerHTML = Tcsi1.toFixed(2);
            tab1.rows[1].cells[1].innerText = (Number(data2.innerHTML)).toFixed(2);
            //温海水出口温度
            document.getElementById("data3").innerHTML = Twso1.toFixed(2);
            tab1.rows[2].cells[1].innerText = (Number(data3.innerHTML)).toFixed(2);
            //冷海水出口温度
            document.getElementById("data4").innerHTML = Tcso1.toFixed(2);
            tab1.rows[3].cells[1].innerText = (Number(data4.innerHTML)).toFixed(2);
            //蒸発温度
            document.getElementById("data5").innerHTML = Te1.toFixed(2);
            tab1.rows[4].cells[1].innerText = (Number(data5.innerHTML)).toFixed(2);
            //凝縮温度
            document.getElementById("data6").innerHTML = Tc1.toFixed(2);
            tab1.rows[5].cells[1].innerText = (Number(data6.innerHTML)).toFixed(2);
            //Point 1 Pressure
            document.getElementById("data7").innerHTML = P[0].toFixed(2);
            tab1.rows[6].cells[1].innerText = (Number(data7.innerHTML)).toFixed(2);
            //Point 2 Pressure
            document.getElementById("data8").innerHTML = P[1].toFixed(2);
            tab1.rows[7].cells[1].innerText = (Number(data8.innerHTML)).toFixed(2);
            //Point 3 Pressure
            document.getElementById("data9").innerHTML = P[2].toFixed(2);
            tab1.rows[8].cells[1].innerText = (Number(data9.innerHTML)).toFixed(2);
            //Point 4 Pressure
            document.getElementById("data10").innerHTML = P[3].toFixed(2);
            tab1.rows[9].cells[1].innerText = (Number(data10.innerHTML)).toFixed(2);
            //Point 1 Temp.
            document.getElementById("data11").innerHTML = T[0].toFixed(2);
            tab1.rows[10].cells[1].innerText = (Number(data11.innerHTML)).toFixed(2);
            //Point 2 Temp.
            document.getElementById("data12").innerHTML = T[1].toFixed(2);
            tab1.rows[11].cells[1].innerText = (Number(data12.innerHTML)).toFixed(2);	
            //Point 3 Temp.
            document.getElementById("data13").innerHTML = T[2].toFixed(2);
            tab1.rows[12].cells[1].innerText = (Number(data13.innerHTML)).toFixed(2);
            //Point 4 Temp.
            document.getElementById("data14").innerHTML = T[3].toFixed(2);
            tab1.rows[13].cells[1].innerText = (Number(data14.innerHTML)).toFixed(2);
            //Turbin output
            document.getElementById("data15").innerHTML = Wt1.toFixed(2);
            tab1.rows[4].cells[4].innerText = (Number(data15.innerHTML)).toFixed(2);
            //冷海水入口温度
            document.getElementById("data17").innerHTML = Tcsi2.toFixed(2);
            tab2.rows[1].cells[1].innerText = (Number(data17.innerHTML)).toFixed(2);
            /*--------低温側--------------*/
            document.getElementById("data16").innerHTML = Twsi2.toFixed(2);
            tab2.rows[0].cells[1].innerText = (Number(data16.innerHTML)).toFixed(2);
            //温海水出口温度
            document.getElementById("data18").innerHTML = Twso2.toFixed(2);
            tab2.rows[2].cells[1].innerText = (Number(data18.innerHTML)).toFixed(2);
            //冷海水出口温度
            document.getElementById("data19").innerHTML = Tcso2.toFixed(2);
            tab2.rows[3].cells[1].innerText = (Number(data19.innerHTML)).toFixed(2);
            //蒸発温度
            document.getElementById("data20").innerHTML = Te2.toFixed(2);
            tab2.rows[4].cells[1].innerText = (Number(data20.innerHTML)).toFixed(2);
            //凝縮温度
            document.getElementById("data21").innerHTML = Tc2.toFixed(2);
            tab2.rows[5].cells[1].innerText = (Number(data21.innerHTML)).toFixed(2);
            //Point 5 Pressure
            document.getElementById("data22").innerHTML = P[4].toFixed(2);
            tab2.rows[6].cells[1].innerText = (Number(data22.innerHTML)).toFixed(2);
            //Point 6 Pressure
            document.getElementById("data23").innerHTML = P[5].toFixed(2);
            tab2.rows[7].cells[1].innerText = (Number(data23.innerHTML)).toFixed(2);
            //Point 7 Pressure
            document.getElementById("data24").innerHTML = P[6].toFixed(2);
            tab2.rows[8].cells[1].innerText = (Number(data24.innerHTML)).toFixed(2);
            //Point 8 Pressure
            document.getElementById("data25").innerHTML = P[7].toFixed(2);
            tab2.rows[9].cells[1].innerText = (Number(data25.innerHTML)).toFixed(2);
            //Point 5 Temp.
            document.getElementById("data26").innerHTML = T[4].toFixed(2);
            tab2.rows[10].cells[1].innerText = (Number(data26.innerHTML)).toFixed(2);
            //Point 6 Temp.
            document.getElementById("data27").innerHTML = T[5].toFixed(2);
            tab2.rows[11].cells[1].innerText = (Number(data27.innerHTML)).toFixed(2);	
            //Point 7 Temp.
            document.getElementById("data28").innerHTML = T[6].toFixed(2);
            tab2.rows[12].cells[1].innerText = (Number(data28.innerHTML)).toFixed(2);
            //Point 8 Temp.
            document.getElementById("data29").innerHTML = T[7].toFixed(2);
            tab2.rows[13].cells[1].innerText = (Number(data29.innerHTML)).toFixed(2);
            //Turbine output
            document.getElementById("data30").innerHTML = Wt2.toFixed(2);
            tab2.rows[4].cells[4].innerText = (Number(data30.innerHTML)).toFixed(2);

            // Data Output of Simulation Results for Graphs of Simulator

            if(haxis == GRAPHDATANUMBER + 1){
                for(i = 1; i <= DATANUMBER; i++){
                    var DataID = "data" + i + ".innerHTML";
                    g[i][GRAPHDATANUMBER + 1] = eval(DataID);
                    g[i].shift();
                }
            }else{
                for(i = 1; i <= DATANUMBER; i++){
                    var DataID = "data" + i + ".innerHTML";
                    g[i][haxis] = eval(DataID);
                }
                haxis = haxis+Number(1);
            }
     
            // Input of Graph Initial value and Step Width
            StartValue=Number(document.js.GStart.value);
            StepWidth=Number(document.js.GStep.value);
            //Time
            presenttime =  (loopNum*DT).toFixed(0)
            document.getElementById("dispT").innerHTML = presenttime;

            tab1.rows[0].cells[4].innerText=(DTme1).toFixed(2);
            tab1.rows[1].cells[4].innerText=(DTmc1).toFixed(2);
            tab1.rows[2].cells[4].innerText=(Twsi1-Twso1).toFixed(2);
            tab1.rows[3].cells[4].innerText=(Tcso1-Tcsi1).toFixed(2);
            tab1.rows[5].cells[4].innerText=(Wp1).toFixed(2);
            tab1.rows[6].cells[4].innerText=(W1).toFixed(2);
            tab1.rows[7].cells[4].innerText=(ETA1).toFixed(2);
            tab1.rows[8].cells[4].innerText=(lastQe1).toFixed(2);
            tab1.rows[9].cells[4].innerText=(lastQc1).toFixed(2);
            tab1.rows[10].cells[4].innerText=(MPR).toFixed(2);

            tab2.rows[0].cells[4].innerText=(DTme2).toFixed(2);
            tab2.rows[1].cells[4].innerText=(DTmc2).toFixed(2);
            tab2.rows[2].cells[4].innerText=(Twsi2-Twso2).toFixed(2);
            tab2.rows[3].cells[4].innerText=(Tcso2-Tcsi2).toFixed(2);
            tab2.rows[5].cells[4].innerText=(Wp2).toFixed(2);
            tab2.rows[6].cells[4].innerText=(W2).toFixed(2);
            tab2.rows[7].cells[4].innerText=(ETA2).toFixed(2);
            tab2.rows[8].cells[4].innerText=(lastQe2).toFixed(2);
            tab2.rows[9].cells[4].innerText=(lastQc2).toFixed(2);
            tab2.rows[10].cells[4].innerText=(MPR).toFixed(2);
        }

        console.log('%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f',(loopNum*DT),Twsi1,Twso1,Te1,Tcsi1,Tcso1,Tc1,Twsi2,Twso2,Te2,Tcsi2,Tcso2,Tc2,Mws,Mcs,Mwf1,Mwf2,DTme1,DTmc1,DTme2,DTmc2,Twsi1-Twso1,Tcso1-Tcsi1,Twsi2-Twso2,Tcso2-Tcsi2,lastQe1,lastQc1,lastQe2,lastQc2,Wt1,Wt2,Wp1,Wp2,W1,W2,ETA1,ETA2,MP,MPR,T[0],P[0],V[0],h[0],S[0],T[1],P[1],V[1],h[1],S[1],T[2],P[2],V[2],h[2],S[2],T[3],P[3],V[3],h[3],S[3],T[4],P[4],V[4],h[4],S[4],T[5],P[5],V[5],h[5],S[5],T[6],P[6],V[6],h[6],S[6],T[7],P[7],V[7],h[7],S[7]);
        timeoutID=setTimeout("mt()",100);    //1秒毎に再起実行
    }

    function sstop(){
        clearTimeout(timeoutID);
        clearTimeout(GtimeoutID);
    }

    function doToggle(UnitID) {
        if(UnitID == 0){
            document.getElementById("bta").style.color="red";
            document.getElementById("btb").style.color="black";
            document.getElementById("btc").style.color="black";

            document.getElementById("TempGraphButtons").style.display="";
            document.getElementById("PressureGraphButtons").style.display="none";
            document.getElementById("FlowGraphButtons").style.display="none";
        }else if(UnitID == 1){
            document.getElementById("bta").style.color="black";
            document.getElementById("btb").style.color="red";
            document.getElementById("btc").style.color="black";

            document.getElementById("TempGraphButtons").style.display="none";
            document.getElementById("PressureGraphButtons").style.display="";
            document.getElementById("FlowGraphButtons").style.display="none";
        }else{
            document.getElementById("bta").style.color="black";
            document.getElementById("btb").style.color="black";
            document.getElementById("btc").style.color="red";

            document.getElementById("TempGraphButtons").style.display="none";
            document.getElementById("PressureGraphButtons").style.display="none";
            document.getElementById("FlowGraphButtons").style.display="";
        }

    }


    GraphTitle = new Array(DATANUMBER + 1);

    GraphTitle[0] = "Graph Title";
    GraphTitle[1] = "Warm Sea Water Inlet T (A) [deg]";
    GraphTitle[2] = "Cold Sea Water Inlet T (A) [deg]";
    GraphTitle[3] = "Warm Sea Water Outlet T (A) [deg]";
    GraphTitle[4] = "Cold Sea Water Outlet T (A) [deg]";
    GraphTitle[5] = "Evaporation Temperature (A) [deg]";
    GraphTitle[6] = "Condensation Temperature (A) [deg]";
    GraphTitle[7] = "Point 1 P [MPa]";
    GraphTitle[8] = "Point 2 P [MPa]";
    GraphTitle[9] = "Point 3 P [MPa]";
    GraphTitle[10] = "Point 4 P [MPa]";
    GraphTitle[11] = "Point 1 T [deg]";
    GraphTitle[12] = "Point 2 T [deg]";
    GraphTitle[13] = "Point 3 T [deg]";
    GraphTitle[14] = "Point 4 T [deg]";
    GraphTitle[15] = "Turbine Output (A) [kW]";
    GraphTitle[16] = "Warm Sea Water Inlet T (B) [deg]";
    GraphTitle[17] = "Cold Sea Water Inlet T (B) [deg]";
    GraphTitle[18] = "Warm Sea Water Outlet T (B) [deg]";
    GraphTitle[19] = "Cold Sea Water Outlet T (B) [deg]";
    GraphTitle[20] = "Evaporation Temperature (B) [deg]";
    GraphTitle[21] = "Condensation Temperature (B) [deg]";
    GraphTitle[22] = "Point 5 P [MPa]";
    GraphTitle[23] = "Point 6 P [MPa]";
    GraphTitle[24] = "Point 7 P [MPa]";
    GraphTitle[25] = "Point 8 P [MPa]";
    GraphTitle[26] = "Point 5 T [deg]";
    GraphTitle[27] = "Point 6 T [deg]";
    GraphTitle[28] = "Point 7 T [deg]";
    GraphTitle[29] = "Point 8 T [deg]";
    GraphTitle[30] = "Turbine Output (B) [kW]";

    //グラフ用関数
    function graphT(VariableID){
        if((1 <= VariableID && VariableID <= 6) || (11 <= VariableID && VariableID <= 14) || (16 <= VariableID && VariableID <= 21) || (26 <= VariableID && VariableID <= 29)){ // Start Value and Step Width for Temperature [DegC]
            StartValue=InitialStartValueTemperature;
            StepWidth=InitialStepWidthTemperature;
            document.getElementById("gstart").value=InitialStartValueTemperature;
            document.getElementById("gstep").value=InitialStepWidthTemperature;
        } else if((7 <= VariableID && VariableID <= 10) || (22 <= VariableID && VariableID <= 25)){ // Start Value and Step Width for Pressure [MPa]
            StartValue=InitialStartValuePressure;
            StepWidth=InitialStepWidthPressure;
            document.getElementById("gstart").value=InitialStartValuePressure;
            document.getElementById("gstep").value=InitialStepWidthPressure;
        } else if(VariableID === 15 || VariableID === 30){ // Start Value and Step Width for Power [kW]
            StartValue=InitialStartValuePower;
            StepWidth=InitialStepWidthPower;
            document.getElementById("gstart").value=InitialStartValuePower;
            document.getElementById("gstep").value=InitialStepWidthPower;
        } else {
            alert("Graph Button Select ERROR!!");
        }
        PreviousVariableID = VariableID;
        document.getElementById("GTitle").value=GraphTitle[VariableID];
    }

    function graph(VariableID){
        if( VariableID == PreviousVariableID){
            var lineChartData = {
                labels : [],
                datasets : [
                    {
                      fillColor : "rgba(0,255,0,0.1)",
                      strokeColor : "rgba(0,255,0,8)",
                      pointColor : "rgba(255,255,255,1)",
                      pointStrokeColor : "#000",
                      data : []
                    }
                ]
            }
            var option = {
                responsive : false,
                scaleOverride : true,
                scaleSteps : 10,
                scaleStepWidth : StepWidth,
                scaleStartValue : StartValue,
                scaleLineColor : "#000000",
                scaleLineWidth : 1,
                scaleShowLabels : true,
                //  scaleFontFamily : "'Arial'",
                scaleFontSize : 10,
                scaleFontStyle : "normal",
                scaleFontColor : "#000000",  
                scaleShowGridLines : false,
                scaleGridLineColor : "#000000",
                scaleGridLineWidth : 1,  
                bezierCurve : false,
                pointDot : true,
                pointDotRadius : 3,
                pointDotStrokeWidth : 1,
                datasetStrokeWidth : 2,
                animation : false,
                animationSteps : 60,
                animationEasing : "easeOutQuart",
            }


            if(presenttime < GRAPHDATANUMBER + 1){
                for(i = 0;i < presenttime;i++){
                    lineChartData.labels[i] = i;
                    lineChartData.datasets[0].data[i] =  g[VariableID][i];
                }
                for(i = presenttime;i < GRAPHDATANUMBER;i++){
                    lineChartData.labels[i] = "";
                    lineChartData.datasets[0].data[i] =  "";
                }
            }else{
                for(i = 0;i < GRAPHDATANUMBER;i++){
                    lineChartData.labels[i] = Number(presenttime) + (i - GRAPHDATANUMBER);
                    lineChartData.datasets[0].data[i] =  g[VariableID][i];
                }
            }


            var ctx = document.getElementById("Gline").getContext("2d");

            ctx.canvas.width = 1250;
            ctx.canvas.height = 370;
            myChart = new Chart(ctx).Line(lineChartData,option);
        }

        var graphname = "graph(" + VariableID + ")";
        GtimeoutID = setTimeout(graphname,1000);
    }