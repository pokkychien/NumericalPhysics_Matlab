%{
    secant_root.m
    割線法找根(假設每個整數之間僅有一個根)
    註: 如果兩個整數之間有超過一個根，那這支沒辦法解。
    Chang Kai-Po @ Jian Lab 2023/03/03
%}

clc; clear; close; 
disp (muller_core(@test_case, 0,1,3))
%fprintf("函數 x^4 + x^3 - 2*x^2 + x - 6 \n 在-20與20之間的根可能有:\n");
%disp(all_root(@test_case, -20, 20, 0.001));

function output = test_case (x)
  %用來測試找根的函數
  output = x^4 + x^3 - 2*x^2 + x - 6;
end 
function output2 = test_case2 (x)
  %用來測試找根的函數
  output2 = sin(3*x)
end 

function output = all_root (f, low, high, epsilon)
   %{
      以low為下界, high為上界，搜尋函數f(x)在這個範圍內可能的根，
      精度為epsilon。      
   %}
    introot = find_integer (f, low, high);      %搜尋函數f(x)在那些整數間可能會有根。
    output = [];
    if isempty(introot)
        return;                                 %如果這範圍沒有，就結束
    end 
    for i = 1:length(introot)
        output(end+1) = one_root (f, introot(i), epsilon);
    end 
end 

function output = one_root (f, i, epsilon)
   %{
      !!除了單元測試用途以外，請不要直接呼叫這個函數!!
      在已經知道函數f在i與i+1之間有一個根的前提下，
      用割線法逼近出這個根在哪裡，精度為epsilon。      
   %}
    low = i; high = i+1;   %一開始的下界是i，上界是i+1
    mid = high - ((high-low)/(f(high)-f(low))*f(high)); %割線公式
    while f(mid)^2 > epsilon        
        %disp([high, mid, low]);
        mid = high - ((high-low)/(f(high)-f(low))*f(high)); %割線公式
        if f(low)*f(mid) < 0            %根在下半，故異號  
            high = mid;                 %中間值變成上界                        
        else                            %根在上半
            low = mid;                  %中間值變成下界
        end 
    end 
    output = mid;                       %以下界為根                    
end 
function output = muller_core (f, x0, x1, x2)
    %{
    !!除了單元測試用途以外，請不要直接呼叫這個函數!!
    給予三個從小到大的x值: x0, x1, x2，並且f(x0), f(x1), f(x2)的符號不同
    算出下一個二次曲線的根
    %}
    h1 = (x1-x0); h2 = (x2-x1);
    delta1 = (f(x1)-f(x0))/h1; delta2 = (f(x2)-f(x1))/h2;
    d = (delta2 - delta1) / (h2 + h1)
    %以上是為了方便運算先設定的一些值
    if delta1 == delta2                            %兩斜率相等，在同一直線上
        output = x2 - ((x2-x1)/(f(x2)-f(x1))*f(x2)); %割線公式
    end 
    b = delta2 + h2*d; 
    D = (b**2 - 4*f(x2)*d)**0.5; 
    if D < 0                                     %不討論虛根
        output = false;
    elif abs(b-D) < abs(b+D)                       %選擇較大的值
        E = b + D;
    else
        E = b - D;
    end 
    h = -2*f(x2)/E;
    output = x2 + h;
 end 

function output = find_integer (f, low, high)
   %{
      以low為下界, high為上界，搜尋函數f(x)在那些整數間可能會有根。
      之後，傳回所有可能會有根的整數下緣。
      (例如說如果有根在2-3之間與7-8之間則傳回[2,7])。
   %}  
    output = [];                        %先放一個空array用來填解
    for i = low:high 
        if f(i)*f(i+1) < 0              %根在這兩個整數之間，故異號
            output(end+1) = i; 
        end 
    end 
end 