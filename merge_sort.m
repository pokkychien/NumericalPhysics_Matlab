clc
clear 
close 

testarray = [1,3,20,2,15,8,9,10];

disp(merge_sort_proper(testarray));

function output = merge_sort_proper(target_array)    
    if length (target_array) ~= 1
        mid = fix (length(target_array)/ 2);
        first_subarray = target_array (1:mid);
        second_subarray = target_array (mid+1:length(target_array));
        first_half = merge_sort_proper (first_subarray);
        second_half = merge_sort_proper (second_subarray);
        output = merge (first_half, second_half);
    else
        output = target_array;
    end 
end 

function combined_array = merge (anterior, posterior)
    p = 1; 
    q = 1;
    combined_array = [];
    while (p <= length(anterior) && q <= length(posterior))
        if anterior(p) <= posterior(q)
            combined_array = [combined_array anterior(p)];
            p = p + 1;                  
        else 
            combined_array = [combined_array posterior(q)];
            q = q + 1;
        end  
    end 
    while (p <= length(anterior))
        combined_array = [combined_array anterior(p)];
        p = p + 1;   
    end 
    while (q <= length(posterior))
        combined_array = [combined_array posterior(q)];
        q = q + 1;   
    end 
   
end 
