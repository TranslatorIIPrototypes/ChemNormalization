kubectl delete deployment cn-neo4j-kp-dep
kubectl delete service cn-neo4j-kp-svc
kubectl delete pvc --all
kubectl delete pv --all
pause
kubectl create -f cn-neo4jkp-kube-conf.yml
kubectl apply -f cn-neo4jkp-persitence.yml
