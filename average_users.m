function avg_users = average_users(pc, lambda_u, pk,pk1,pk2,pk3,l1,l2,l3, y)

avg_users = (pc*lambda_u*power(pk,2/y))/(l1*power(pk1,2/y)+l2*power(pk2,2/y)+l3*power(pk3,2/y));

end