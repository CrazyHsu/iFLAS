$(document).ready(function(){
    var d=0;
    var a=false;
    var c=new Array();
    var b=new Array();
    $(document).keyup(function(g){
      var f=new Date().getTime();
      if(g.keyCode===17){
        var h=f-d;
        d=f;
        if(h<500){
          if(a){
            $(".cb-search-tool").css("display","none");
            a=false}
          else{
            $(".cb-search-tool").css("display","block");
            a=true;
            $("#cb-search-content").val("");
            $("#cb-search-content").focus()}
          d=0}}
      else{
        if(g.keyCode==27){
          $(".cb-search-tool").css("display","none");
          a=false;
          d=0}
      }
    });
    $("#cb-search-content").keyup(function(g){
      var f=new Date().getTime();
      if(g.keyCode==17){
        var h=f-d;
        d=f;
        if(h<500){
          if(a){
            $(".cb-search-tool").css("display","none");
            a=false
          }else{
            $(".cb-search-tool").css("display","block");
            a=true;
            $("#cb-search-content").val("");
            $("#cb-search-content").focus()
          }
          d=0
        }
      }
    });
    $("#cb-close-btn").click(function(){
      $(".cb-search-tool").css("display","none");
      a=false;
      d=0
    });
    $("#cb-search-btn").click(function(){
      $(".cb-search-tool").css("display","block");
      a=true;
      $("#cb-search-content").val("");
      $("#cb-search-content").focus();
      d=0
    });
    $.getJSON("/search/cb-search.json").done(function(g){
      if(g.code==0){
        for(var e in g.data){
          var f=g.data[e];
          c.push(f.title);
          b.push(f.url)
        }
        $("#cb-search-content").typeahead({source:c,afterSelect:function(h){
          $(".cb-search-tool").css("display","none");
          a=false;
          window.location.href=(b[c.indexOf(h)]);
          return h}
        })
      }
    })
});