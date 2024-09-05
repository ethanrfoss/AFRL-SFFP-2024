
function v = QuaternionRotate(q,v)

v = [zeros(3,1) eye(3)]*QuaternionMultiply(QuaternionMultiply(q,[0;v]),QuaternionConjugate(q));

end