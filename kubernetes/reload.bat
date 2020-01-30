kubectl delete -f .
pause

kubectl create -f cn-neo4jkp-kube-conf.yml
kubectl apply -f cn-neo4jkp-persitence.yml
minikube service cn-neo4j-kp-svc --url
