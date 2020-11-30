let params = {
    d:9,
    i:15,
    s:3,
    c:15
}

function maxKeys(obj){
let a = [];

for (var i in obj) {
let v = obj[i];
a[v] || (a[v] = []);
a[v].push(i);
}

return a.pop();
}


arr = maxKeys(params)
console.log(arr)

