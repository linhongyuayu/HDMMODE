function outcome = Convertion2(Average,G,last_gen,change_threshold)
    if ( max(abs(Average(G,:) - Average(G-last_gen,:)) ) <= change_threshold)
        outcome = true;
    else
        outcome = false;
    end
end