$(document).ready(main);
var contador = 1;
function main(){
    $('.menu').click(function(){
        if(contador==1){
        $('.b').animate({
            left: '+41%'
        });
        contador=0;
    }else{
        contador=1;
        $('.b').animate({
            left: '-41%'
        });
    }
})
}